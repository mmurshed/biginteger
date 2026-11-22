#pragma once

namespace BigMath::detail
{
  // Multiplication can be reached through several layers of Newton division.
  // Each conversion subtree sets the policy on the thread doing its arithmetic.
  class DecimalConversionScope
  {
  public:
    DecimalConversionScope() noexcept : previous(active) { active = true; }
    ~DecimalConversionScope() { active = previous; }
    DecimalConversionScope(const DecimalConversionScope &) = delete;
    DecimalConversionScope &operator=(const DecimalConversionScope &) = delete;

    static bool Active() noexcept { return active; }

  private:
    static inline thread_local bool active = false;
    bool previous;
  };
}
