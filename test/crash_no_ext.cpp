#include <cstdint>

#include <fbitset.hpp>

int main()
{
    fbitset::Fbitset<1, uint32_t, fbitset::No_ext> bits(33);
    return 0;
}
