#pragma once

#include "vcl/vcl.hpp"
#include <random>
#include <functional>

using namespace vcl;

class PerlinNoise
{
    static const unsigned tableSize = 256;
    static const unsigned tableSizeMask = tableSize - 1;
    vec3 gradients[tableSize];
    unsigned permutationTable[tableSize * 2];


    inline
    float lerp(const float lo, const float hi, const float t) const {
        return lo * (1 - t) + hi * t;
    }

    inline
    float smoothstep(const float &t) const {
        return t * t * (3 - 2 * t);
    }

    /* inline */
    int hash(const int &x, const int &y, const int &z) const
    { return permutationTable[permutationTable[permutationTable[x] + y] + z]; }

public:
    PerlinNoise()
    {
        unsigned seed = 2016;
        std::mt19937 generator(seed);
        std::uniform_real_distribution<float> distribution;
        auto dice = std::bind(distribution, generator);

        for (unsigned i = 0; i < tableSize; ++i) {
            float theta = acos(2 * dice() - 1);
            float phi = 2 * dice() * M_PI;

            float x = cos(phi) * sin(theta);
            float y = sin(phi) * sin(theta);
            float z = cos(theta);
            gradients[i] = vec3(x, y, z);
        }

        std::uniform_int_distribution<unsigned> distributionInt;
        auto diceInt = std::bind(distributionInt, generator);
        // create permutation table
        for (unsigned i = 0; i < tableSize; ++i)
            std::swap(permutationTable[i], permutationTable[diceInt() & tableSizeMask]);
        // extend the permutation table in the index range [256:512]
        for (unsigned i = 0; i < tableSize; ++i) {
            permutationTable[i] &= tableSizeMask;
            permutationTable[tableSize + i] = permutationTable[i];
        }
    }
    virtual ~PerlinNoise() {}

    float eval(const vec3 &p) const
    {
        int xi0 = ((int)std::floor(p.x)) & tableSizeMask;
        int yi0 = ((int)std::floor(p.y)) & tableSizeMask;
        int zi0 = ((int)std::floor(p.z)) & tableSizeMask;

        int xi1 = (xi0 + 1) & tableSizeMask;
        int yi1 = (yi0 + 1) & tableSizeMask;
        int zi1 = (zi0 + 1) & tableSizeMask;

        float tx = p.x - ((int)std::floor(p.x));
        float ty = p.y - ((int)std::floor(p.y));
        float tz = p.z - ((int)std::floor(p.z));

        float u = smoothstep(tx);
        float v = smoothstep(ty);
        float w = smoothstep(tz);

        // gradients at the corner of the cell
        const vec3 &c000 = gradients[hash(xi0, yi0, zi0)];
        const vec3 &c100 = gradients[hash(xi1, yi0, zi0)];
        const vec3 &c010 = gradients[hash(xi0, yi1, zi0)];
        const vec3 &c110 = gradients[hash(xi1, yi1, zi0)];

        const vec3 &c001 = gradients[hash(xi0, yi0, zi1)];
        const vec3 &c101 = gradients[hash(xi1, yi0, zi1)];
        const vec3 &c011 = gradients[hash(xi0, yi1, zi1)];
        const vec3 &c111 = gradients[hash(xi1, yi1, zi1)];

        // generate vectors going from the grid points to p
        float x0 = tx, x1 = tx - 1;
        float y0 = ty, y1 = ty - 1;
        float z0 = tz, z1 = tz - 1;

        vec3 p000 = vec3(x0, y0, z0);
        vec3 p100 = vec3(x1, y0, z0);
        vec3 p010 = vec3(x0, y1, z0);
        vec3 p110 = vec3(x1, y1, z0);

        vec3 p001 = vec3(x0, y0, z1);
        vec3 p101 = vec3(x1, y0, z1);
        vec3 p011 = vec3(x0, y1, z1);
        vec3 p111 = vec3(x1, y1, z1);

        // linear interpolation
        float a = lerp(dot(c000, p000), dot(c100, p100), u);
        float b = lerp(dot(c010, p010), dot(c110, p110), u);
        float c = lerp(dot(c001, p001), dot(c101, p101), u);
        float d = lerp(dot(c011, p011), dot(c111, p111), u);

        float e = lerp(a, b, v);
        float f = lerp(c, d, v);

        return lerp(e, f, w); // g
    }
};
