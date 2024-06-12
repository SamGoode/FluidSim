#include "RadixSort.h"

// radix sort is O(d(n + w)) where w is size of digit and d is number of digits
// for 128 different cell IDs I want to find the sweet spot
// w = 16, (base 16) 4 bits
// 127 = 0x7F
// 2 digits
// w = 10 (decimal)
// 3 digits
// w = 2, (binary) 1 bit
// 127 = 01111111
// 7 digits

// radix sort
void radixSort(Array<int>& unsorted) {
    int maxValue = 0;
    for (int i = 0; i < unsorted.getCount(); i++) {
        if (maxValue < unsorted[i]) {
            maxValue = unsorted[i];
        }
    }

    Array<int> sorted(unsorted.getCount(), unsorted.getCapacity());
    for (int i = 0; (1 << (i * 4)) < maxValue; i++) {
        // contains the range of possible values found in a digit
        Array<int> count(16, 16, 0);

        for (int n = 0; n < unsorted.getCount(); n++) {
            char halfByte = (unsorted[n] >> (i * 4)) & 0x0000000F;
            count[halfByte]++;
        }

        for (int n = 1; n < count.getCount(); n++) {
            count[n] += count[n - 1];
        }

        for (int n = unsorted.getCount() - 1; n >= 0; n--) {
            char halfByte = (unsorted[n] >> (i * 4)) & 0x0000000F;
            sorted[count[halfByte] - 1] = unsorted[n];
            count[halfByte]--;
        }

        unsorted.quickCopy(sorted);
    }
}

void radixSort(Array<int2>& unsorted) {
    int maxValue = 0;
    for (int i = 0; i < unsorted.getCount(); i++) {
        if (maxValue < unsorted[i].y) {
            maxValue = unsorted[i].y;
        }
    }

    Array<int2> sorted(unsorted.getCount(), unsorted.getCapacity());
    for (int i = 0; (1 << (i * 4)) < maxValue; i++) {
        // contains the range of possible values found in a digit
        Array<int> count(16, 16, 0);

        for (int n = 0; n < unsorted.getCount(); n++) {
            char halfByte = (unsorted[n].y >> (i * 4)) & 0x0000000F;
            count[halfByte]++;
        }

        for (int n = 1; n < count.getCount(); n++) {
            count[n] += count[n - 1];
        }

        for (int n = unsorted.getCount() - 1; n >= 0; n--) {
            char halfByte = (unsorted[n].y >> (i * 4)) & 0x0000000F;
            sorted[count[halfByte] - 1] = unsorted[n];
            count[halfByte]--;
        }

        unsorted.quickCopy(sorted);
    }
}

void radixSort(Array<int2>& unsorted, int maxValue) {
    Array<int2> sorted(unsorted.getCount(), unsorted.getCapacity());
    for (int i = 0; (1 << (i * 4)) < maxValue; i++) {
        // contains the range of possible values found in a digit
        Array<int> count(16, 16, 0);

        for (int n = 0; n < unsorted.getCount(); n++) {
            char halfByte = (char)(unsorted[n].y >> (i * 4)) & 0x0F;
            count[halfByte]++;
        }

        for (int n = 1; n < count.getCount(); n++) {
            count[n] += count[n - 1];
        }

        for (int n = unsorted.getCount() - 1; n >= 0; n--) {
            char halfByte = (char)(unsorted[n].y >> (i * 4)) & 0x0F;
            sorted[count[halfByte] - 1] = unsorted[n];
            count[halfByte]--;
        }

        unsorted.quickCopy(sorted);
    }
}


// count sort
Array<int> countSort(Array<int> unsorted) {
    // contains the range of possible values found in unsorted array
    Array<int> count(16, 0);
    for (int i = 0; i < unsorted.getCount(); i++) {
        count[unsorted[i]] += 1;
    }

    for (int i = 1; i < count.getCount(); i++) {
        count[i] += count[i - 1];
    }

    Array<int> sorted(unsorted.getCount());
    for (int i = unsorted.getCount() - 1; i >= 0; i--) {
        sorted[count[unsorted[i]] - 1] = unsorted[i];
        count[unsorted[i]] -= 1;
    }

    return sorted;
}