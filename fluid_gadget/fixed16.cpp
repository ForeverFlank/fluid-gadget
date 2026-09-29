#include <inttypes.h>

const int16_t SQRT_LOOKUP[] = {
    0, 256, 362, 443, 512, 572, 627, 677,
    724, 768, 809, 849, 886, 923, 957, 991,
    1024, 1055, 1086, 1115, 1144, 1173, 1200, 1227,
    1254, 1280, 1305, 1330, 1354, 1378, 1402, 1425,
    1448, 1470, 1492, 1514, 1536, 1557, 1578, 1598,
    1619, 1639, 1659, 1678, 1698, 1717, 1736, 1755,
    1773, 1792, 1810, 1828, 1846, 1863, 1881, 1898,
    1915, 1932, 1949, 1966, 1982, 1999, 2015, 2031,
    2048, 2063, 2079, 2095, 2111, 2126, 2141, 2157,
    2172, 2187, 2202, 2217, 2231, 2246, 2260, 2275,
    2289, 2304, 2318, 2332, 2346, 2360, 2374, 2387,
    2401, 2415, 2428, 2442, 2455, 2468, 2482, 2495,
    2508, 2521, 2534, 2547, 2560, 2572, 2585, 2598,
    2610, 2623, 2635, 2648, 2660, 2672, 2684, 2697,
    2709, 2721, 2733, 2745, 2757, 2769, 2780, 2792,
    2804, 2816, 2827, 2839, 2850, 2862, 2873, 2884,
};

struct fixed16
{
    int16_t value;


    fixed16()
    {
        this->value = 0;
    }

    fixed16(float value)
    {
        this->value = (int16_t)(value * 256.f);
    }

    fixed16(int16_t value)
    {
        this->value = value;
    }


    float to_float()
    {
        return value / 256.f;
    }

    int to_int()
    {
        return value >> 8;
    }


    fixed16 operator+()
    {
        return fixed16(static_cast<int16_t>(+value));
    }

    fixed16 operator-()
    {
        return fixed16(static_cast<int16_t>(-value));
    }

    fixed16 operator+(fixed16 rhs)
    {
        return fixed16(static_cast<int16_t>(value + rhs.value));
    }

    fixed16 operator-(fixed16 rhs)
    {
        return fixed16(static_cast<int16_t>(value - rhs.value));
    }

    fixed16 operator*(fixed16 rhs)
    {
        int32_t res = value * rhs.value;
        return fixed16(static_cast<int16_t>(res >> 8));
    }

    fixed16 operator/(fixed16 rhs)
    {
        int32_t lhs32 = value << 8;
        return fixed16(static_cast<int16_t>(lhs32 / rhs.value));
    }


    void operator+=(fixed16 rhs) { (*this) = (*this) + rhs; }

    void operator-=(fixed16 rhs) { (*this) = (*this) - rhs; }

    void operator*=(fixed16 rhs) { (*this) = (*this) * rhs; }

    void operator/=(fixed16 rhs) { (*this) = (*this) / rhs; }


    bool operator<(fixed16 rhs) { return value < rhs.value; }

    bool operator>(fixed16 rhs) { return value > rhs.value; }


    // operator float() const { return  this->value / 256.f; }
};


fixed16 sqrt(fixed16 value)
{
    if (value.value < 0) return fixed16(0.f);

    int8_t upper = value.value >> 8;
    int8_t lower = value.value & 0xFF;

    fixed16 l = fixed16(static_cast<int16_t>(SQRT_LOOKUP[upper]));
    fixed16 h = fixed16(static_cast<int16_t>(SQRT_LOOKUP[upper + 1]));

    fixed16 c = fixed16(static_cast<int16_t>(lower)) / fixed16(1.f);

    return l * (fixed16(1.f) - c) + h * c;
}

fixed16 floor(fixed16 value)
{
    return fixed16(static_cast<int16_t>(value.value & 0xFF00));
}

fixed16 min(fixed16 lhs, fixed16 rhs)
{
    return lhs < rhs ? lhs : rhs;
}

fixed16 max(fixed16 lhs, fixed16 rhs)
{
    return lhs > rhs ? lhs : rhs;
}