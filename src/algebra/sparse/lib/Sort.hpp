#pragma once
#include "algebra/core/types.hpp"

namespace karu {

template<typename index_type>
index_type getMax(index_type keys[], u32 n)
{
    index_type mx = keys[0];
    for (u32 i = 1; i < n; i++)
        if (keys[i] > mx)
            mx = keys[i];
    return mx;
}

template<typename index_type, typename value_type>
void countSort(index_type keys[], value_type vals[], u32 n, u32 exp)
{
    index_type out_keys[n];
    value_type out_vals[n];
	
    i32 i, count[10] = { 0 };
 
    for (i = 0; i < n; i++)
        count[(keys[i] / exp) % 10]++;
 
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    for (i = n - 1; i >= 0; i--) {
        out_keys[count[(keys[i] / exp) % 10] - 1] = keys[i];
        out_vals[count[(keys[i] / exp) % 10] - 1] = vals[i];
        count[(keys[i] / exp) % 10]--;
    }
 
    for (i = 0; i < n; i++)
		{
			keys[i] = out_keys[i];
			vals[i] = out_vals[i];
		}
}

template<typename index_type, typename value_type>
void sortKeyValuePair(index_type keys[], value_type vals[], u32 n)
{

    u32 m = getMax(keys, n);

    for (i32 exp = 1; m / exp > 0; exp *= 10)
        countSort(keys, vals, n, exp);
}

}

 