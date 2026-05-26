/**
 * BigMath calculator REPL.
 *
 * Production interactive frontend over ExpressionEvaluator + BigInteger.
 * Features:
 *   - prompt only when stdin is a TTY (clean for piping)
 *   - persistent variable bindings; implicit `ans` for the last result
 *   - directives: :help :quit :vars :reset :base :digits :time :load :save
 *   - decimal, hexadecimal, and binary output bases
 *   - multi-line continuation via trailing '\\'
 *   - '#' line comments
 *   - structured errors with caret pointing at the offending byte
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/IO.h"
#include "ExpressionEvaluator.h"

using namespace BigMath;

namespace
{
  constexpr char const *kVersion = "12.0";

  struct Settings
  {
    int  base       = 10;     // 10, 16, or 2
    int  maxDigits  = 0;      // 0 = unlimited
    bool showTime   = false;
    bool isTty      = false;
  };

  // ─── output formatting ──────────────────────────────────────────────────────

  std::string FormatHex(BigInteger const &x)
  {
    std::string out;
    if (x.IsNegative()) out.push_back('-');
    out += "0x";
    auto const &v = x.GetInteger();
    int hi = (int)v.size() - 1;
    while (hi > 0 && v[hi] == 0) --hi;
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%x", (unsigned)v[hi]);
    out += buf;
    for (int i = hi - 1; i >= 0; --i)
    {
      std::snprintf(buf, sizeof(buf), "%08x", (unsigned)v[i]);
      out += buf;
    }
    return out;
  }

  std::string FormatBin(BigInteger const &x)
  {
    std::string out;
    if (x.IsNegative()) out.push_back('-');
    out += "0b";
    auto const &v = x.GetInteger();
    int hi = (int)v.size() - 1;
    while (hi > 0 && v[hi] == 0) --hi;

    auto pushBits = [&](uint32_t w, int bits) {
      for (int i = bits - 1; i >= 0; --i)
        out.push_back(((w >> i) & 1u) ? '1' : '0');
    };

    // Top limb: strip leading zero bits (but keep at least one).
    uint32_t top = (uint32_t)v[hi];
    int topBits = 32;
    while (topBits > 1 && ((top >> (topBits - 1)) & 1u) == 0)
      --topBits;
    pushBits(top, topBits);
    for (int i = hi - 1; i >= 0; --i)
      pushBits((uint32_t)v[i], 32);
    return out;
  }

  // Truncate the middle of `s` to fit within maxDigits characters of digits
  // (sign and base prefix are not counted). 0 means no truncation.
  std::string MaybeTruncate(std::string s, int maxDigits)
  {
    if (maxDigits <= 0) return s;
    size_t prefix = 0;
    if (!s.empty() && s[0] == '-') prefix = 1;
    if (s.size() > prefix + 2 && s[prefix] == '0' && (s[prefix + 1] == 'x' || s[prefix + 1] == 'b'))
      prefix += 2;
    size_t body = s.size() - prefix;
    if ((int)body <= maxDigits) return s;
    int head = maxDigits / 2;
    int tail = maxDigits - head;
    std::string out;
    out.reserve(prefix + maxDigits + 16);
    out.append(s, 0, prefix + head);
    out.append("…(");
    out += std::to_string(body - head - tail);
    out.append(" more)…");
    out.append(s, s.size() - tail, tail);
    return out;
  }

  std::string Format(BigInteger const &x, Settings const &s)
  {
    std::string raw;
    switch (s.base)
    {
      case 16: raw = FormatHex(x); break;
      case 2:  raw = FormatBin(x); break;
      default: raw = ToString(x);  break;
    }
    return MaybeTruncate(std::move(raw), s.maxDigits);
  }

  // ─── small string helpers ───────────────────────────────────────────────────

  std::string Trim(std::string const &s)
  {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
  }

  bool IsIdentChar(char c, bool first)
  {
    if (c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) return true;
    return !first && (c >= '0' && c <= '9');
  }

  // Detect `name = expr` at top level (no '=' inside parens).
  // Returns true and fills `name` / `expr` if matched. The `==` of equality
  // isn't supported by the evaluator so we don't need to disambiguate it.
  bool DetectAssignment(std::string const &line, std::string &name, std::string &expr)
  {
    int depth = 0;
    for (size_t i = 0; i < line.size(); ++i)
    {
      char c = line[i];
      if (c == '(') ++depth;
      else if (c == ')') --depth;
      else if (c == '=' && depth == 0)
      {
        std::string left = Trim(line.substr(0, i));
        if (left.empty()) return false;
        for (size_t j = 0; j < left.size(); ++j)
          if (!IsIdentChar(left[j], j == 0)) return false;
        name = left;
        expr = Trim(line.substr(i + 1));
        return true;
      }
    }
    return false;
  }

  std::string CaretLine(std::string const &src, int pos)
  {
    std::string out;
    out.reserve(src.size() + (size_t)pos + 4);
    out += "    ";
    out += src;
    out += "\n    ";
    out.append((size_t)std::max(0, pos), ' ');
    out += '^';
    return out;
  }

  // ─── REPL loop ──────────────────────────────────────────────────────────────

  void PrintHelp()
  {
    std::cout <<
        "BigMath calculator " << kVersion << "\n"
        "\n"
        "Expressions:\n"
        "  + - * / %                arithmetic (integer)\n"
        "  ^                        exponent (right-assoc, non-negative)\n"
        "  ( )                      grouping\n"
        "  unary + -                e.g. -(-5) == 5\n"
        "  123_456_789              underscore digit separators\n"
        "  ans                      last computed value\n"
        "  name = expr              assign to a variable\n"
        "  # comment to end of line\n"
        "  trailing backslash \\     continue on the next line\n"
        "\n"
        "Directives:\n"
        "  :help              show this help\n"
        "  :quit              exit (also Ctrl-D / EOF)\n"
        "  :vars              list variables\n"
        "  :reset             clear all variables\n"
        "  :base dec|hex|bin  output base (default: dec)\n"
        "  :digits N          truncate output to N digits (0 = unlimited)\n"
        "  :time on|off       show evaluation time\n"
        "  :load NAME PATH    read variable from file (decimal text)\n"
        "  :save NAME PATH    write variable to file (decimal text)\n"
        ;
  }

  void PrintBanner()
  {
    std::cout << "BigMath calculator " << kVersion
              << "  (type :help for commands, :quit or Ctrl-D to exit)\n";
  }

  void ListVars(std::map<std::string, BigInteger> const &vars, Settings const &s)
  {
    if (vars.empty())
    {
      std::cout << "(no variables)\n";
      return;
    }
    for (auto const &[k, v] : vars)
      std::cout << "  " << k << " = " << Format(v, s) << "\n";
  }

  bool RunDirective(std::string const &line,
                    Settings &s,
                    std::map<std::string, BigInteger> &vars,
                    bool &shouldQuit)
  {
    // tokenise by whitespace
    std::vector<std::string> tok;
    std::istringstream is(line);
    std::string t;
    while (is >> t) tok.push_back(t);
    if (tok.empty()) return true;

    std::string const &cmd = tok[0];
    if (cmd == ":help" || cmd == ":h" || cmd == ":?")    { PrintHelp(); return true; }
    if (cmd == ":quit" || cmd == ":q" || cmd == ":exit") { shouldQuit = true; return true; }
    if (cmd == ":vars")                                   { ListVars(vars, s); return true; }
    if (cmd == ":reset")
    {
      vars.clear();
      std::cout << "Cleared variables.\n";
      return true;
    }
    if (cmd == ":base")
    {
      if (tok.size() != 2) { std::cerr << "usage: :base dec|hex|bin\n"; return true; }
      if      (tok[1] == "dec" || tok[1] == "10") s.base = 10;
      else if (tok[1] == "hex" || tok[1] == "16") s.base = 16;
      else if (tok[1] == "bin" || tok[1] == "2")  s.base = 2;
      else { std::cerr << "unknown base: " << tok[1] << "\n"; return true; }
      std::cout << "base = " << s.base << "\n";
      return true;
    }
    if (cmd == ":digits")
    {
      if (tok.size() != 2) { std::cerr << "usage: :digits N\n"; return true; }
      try { s.maxDigits = std::stoi(tok[1]); }
      catch (...) { std::cerr << "bad N\n"; return true; }
      if (s.maxDigits < 0) s.maxDigits = 0;
      std::cout << "digits = " << (s.maxDigits ? std::to_string(s.maxDigits) : std::string("unlimited")) << "\n";
      return true;
    }
    if (cmd == ":time")
    {
      if (tok.size() != 2) { std::cerr << "usage: :time on|off\n"; return true; }
      if      (tok[1] == "on")  s.showTime = true;
      else if (tok[1] == "off") s.showTime = false;
      else { std::cerr << "expected on|off\n"; return true; }
      std::cout << "time = " << (s.showTime ? "on" : "off") << "\n";
      return true;
    }
    if (cmd == ":load")
    {
      if (tok.size() != 3) { std::cerr << "usage: :load NAME PATH\n"; return true; }
      std::ifstream in(tok[2]);
      if (!in) { std::cerr << "cannot open: " << tok[2] << "\n"; return true; }
      std::string body((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
      try { vars[tok[1]] = BigIntegerBuilder::From(Trim(body)); }
      catch (...) { std::cerr << "parse failed\n"; return true; }
      std::cout << tok[1] << " loaded (" << vars[tok[1]].size() << " limbs)\n";
      return true;
    }
    if (cmd == ":save")
    {
      if (tok.size() != 3) { std::cerr << "usage: :save NAME PATH\n"; return true; }
      auto it = vars.find(tok[1]);
      if (it == vars.end()) { std::cerr << "no such variable: " << tok[1] << "\n"; return true; }
      std::ofstream out(tok[2]);
      if (!out) { std::cerr << "cannot write: " << tok[2] << "\n"; return true; }
      out << ToString(it->second);
      std::cout << "wrote " << tok[2] << "\n";
      return true;
    }

    std::cerr << "unknown directive: " << cmd << " (try :help)\n";
    return true;
  }

  bool ReadLineWithContinuation(std::istream &in, std::string &out)
  {
    if (!std::getline(in, out)) return false;
    while (!out.empty() && out.back() == '\\')
    {
      out.pop_back();
      std::string cont;
      if (!std::getline(in, cont)) break;
      out += cont;
    }
    return true;
  }
}

int main(int argc, char **argv)
{
  Settings settings;
  settings.isTty = isatty(fileno(stdin));

  for (int i = 1; i < argc; ++i)
  {
    std::string a = argv[i];
    if (a == "--help" || a == "-h") { PrintHelp(); return 0; }
    if (a == "--version" || a == "-V") { std::cout << kVersion << "\n"; return 0; }
    std::cerr << "unknown option: " << a << "\n";
    return 2;
  }

  std::map<std::string, BigInteger> vars;
  ExpressionEvaluator::Lookup lookup =
      [&vars](std::string const &name, BigInteger &out) -> bool
      {
        auto it = vars.find(name);
        if (it == vars.end()) return false;
        out = it->second;
        return true;
      };
  ExpressionEvaluator eval(lookup);

  if (settings.isTty) PrintBanner();

  std::string line;
  bool shouldQuit = false;
  while (!shouldQuit)
  {
    if (settings.isTty) std::cout << "> " << std::flush;
    if (!ReadLineWithContinuation(std::cin, line))
    {
      if (settings.isTty) std::cout << "\n";
      break;
    }

    std::string trimmed = Trim(line);
    if (trimmed.empty() || trimmed[0] == '#') continue;

    if (trimmed[0] == ':')
    {
      RunDirective(trimmed, settings, vars, shouldQuit);
      continue;
    }

    std::string varName, exprText;
    bool assigning = DetectAssignment(trimmed, varName, exprText);
    if (!assigning) exprText = trimmed;

    if (exprText.empty())
    {
      std::cerr << "empty expression\n";
      continue;
    }

    try
    {
      auto t0 = std::chrono::steady_clock::now();
      BigInteger result = eval.Eval(exprText);
      auto t1 = std::chrono::steady_clock::now();

      vars["ans"] = result;
      if (assigning) vars[varName] = result;

      std::cout << Format(result, settings) << "\n";

      if (settings.showTime)
      {
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "(" << ms << " ms)\n";
      }
    }
    catch (EvalError const &e)
    {
      std::cerr << "error: " << e.what() << "\n"
                << CaretLine(exprText, e.Position()) << "\n";
    }
    catch (std::exception const &e)
    {
      std::cerr << "error: " << e.what() << "\n";
    }
    catch (char const *msg)
    {
      std::cerr << "error: " << msg << "\n";
    }
    catch (...)
    {
      std::cerr << "error: unknown exception\n";
    }
  }

  return 0;
}
