#pragma once
// math_compat.h
// Make long-double math functions visible in namespace std if the platform header
// didn't already put them there. This header is intended to be force-included
// via the compiler -include option.

#include <cmath>

// If the C declarations are not visible in global namespace, declare them here.
// Use extern "C" to match C linkage for the math library.
extern "C" {
    long double expl(long double) noexcept;
    long double logl(long double) noexcept;
    long double sqrtl(long double) noexcept;
    long double log1pl(long double) noexcept;
    long double expm1l(long double) noexcept;
}

// Pull these into std:: so code that uses std::expl/std::logl compiles.
namespace std {
    using ::expl;
    using ::logl;
    using ::sqrtl;
    using ::log1pl;
    using ::expm1l;
}
