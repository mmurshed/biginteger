// Minimal in-tree unit-test framework for BigMath.
// No external dependencies; one binary aggregates all REGISTER_TEST entries
// from every translation unit via a static-initializer registry.

#ifndef BIGMATH_UNIT_TEST_FRAMEWORK_H
#define BIGMATH_UNIT_TEST_FRAMEWORK_H

#include <concepts>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace bigmath_ut
{
  struct TestCase
  {
    const char *suite;
    const char *name;
    std::function<void()> fn;
  };

  inline std::vector<TestCase> &Registry()
  {
    static std::vector<TestCase> r;
    return r;
  }

  struct Registrar
  {
    Registrar(const char *suite, const char *name, std::function<void()> fn)
    {
      Registry().push_back({suite, name, std::move(fn)});
    }
  };

  struct AssertionFailure
  {
    std::string message;
  };

  // Concept: T is streamable if std::ostream << T compiles.
  template <typename T>
  concept Streamable = requires(std::ostream &os, const T &v) {
    { os << v } -> std::same_as<std::ostream &>;
  };

  template <Streamable T>
  inline std::string ToStr(const T &v)
  {
    std::ostringstream os;
    os << v;
    return os.str();
  }

  template <typename T>
  inline std::string ToStr(const T &)
  {
    return "<non-streamable>";
  }

  inline void FailAt(const char *file, int line, const std::string &msg)
  {
    std::ostringstream os;
    os << file << ":" << line << ": " << msg;
    throw AssertionFailure{os.str()};
  }
}

// ── assertion macros ─────────────────────────────────────────────────────────

#define ASSERT_TRUE(expr)                                                                              \
  do {                                                                                                 \
    if (!(expr))                                                                                       \
      ::bigmath_ut::FailAt(__FILE__, __LINE__, "ASSERT_TRUE(" #expr ") failed");                       \
  } while (0)

#define ASSERT_FALSE(expr)                                                                             \
  do {                                                                                                 \
    if ((expr))                                                                                        \
      ::bigmath_ut::FailAt(__FILE__, __LINE__, "ASSERT_FALSE(" #expr ") failed");                      \
  } while (0)

#define ASSERT_EQ(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l == _r))                                                                                   \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_EQ(" #lhs ", " #rhs ") failed: " + ::bigmath_ut::ToStr(_l) +                         \
          " != " + ::bigmath_ut::ToStr(_r));                                                           \
  } while (0)

#define ASSERT_NE(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l != _r))                                                                                   \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_NE(" #lhs ", " #rhs ") failed: both = " + ::bigmath_ut::ToStr(_l));                  \
  } while (0)

#define ASSERT_LT(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l < _r))                                                                                    \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_LT(" #lhs ", " #rhs ") failed: " + ::bigmath_ut::ToStr(_l) +                         \
          " not < " + ::bigmath_ut::ToStr(_r));                                                        \
  } while (0)

#define ASSERT_LE(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l <= _r))                                                                                   \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_LE(" #lhs ", " #rhs ") failed: " + ::bigmath_ut::ToStr(_l) +                         \
          " not <= " + ::bigmath_ut::ToStr(_r));                                                       \
  } while (0)

#define ASSERT_GT(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l > _r))                                                                                    \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_GT(" #lhs ", " #rhs ") failed: " + ::bigmath_ut::ToStr(_l) +                         \
          " not > " + ::bigmath_ut::ToStr(_r));                                                        \
  } while (0)

#define ASSERT_GE(lhs, rhs)                                                                            \
  do {                                                                                                 \
    auto &&_l = (lhs);                                                                                 \
    auto &&_r = (rhs);                                                                                 \
    if (!(_l >= _r))                                                                                   \
      ::bigmath_ut::FailAt(__FILE__, __LINE__,                                                         \
          "ASSERT_GE(" #lhs ", " #rhs ") failed: " + ::bigmath_ut::ToStr(_l) +                         \
          " not >= " + ::bigmath_ut::ToStr(_r));                                                       \
  } while (0)

// Token-pasting helpers
#define BIGMATH_UT_CAT2(a, b) a##b
#define BIGMATH_UT_CAT(a, b)  BIGMATH_UT_CAT2(a, b)

// REGISTER_TEST(suite, name) { body }
#define REGISTER_TEST(suite, name)                                                                     \
  static void BIGMATH_UT_CAT(suite, BIGMATH_UT_CAT(_, name))();                                        \
  static ::bigmath_ut::Registrar BIGMATH_UT_CAT(reg_, BIGMATH_UT_CAT(suite, BIGMATH_UT_CAT(_, name)))( \
      #suite, #name, &BIGMATH_UT_CAT(suite, BIGMATH_UT_CAT(_, name)));                                 \
  static void BIGMATH_UT_CAT(suite, BIGMATH_UT_CAT(_, name))()

#endif
