#pragma once
#include "Array.h"
#include "int2.h"

// radix sort
void radixSort(Array<int>& unsorted);
void radixSort(Array<int2>& unsorted);
void radixSort(Array<int2>& unsorted, int maxValue);

// count sort
Array<int> countSort(Array<int> unsorted);