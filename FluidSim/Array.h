#pragma once
#include <string>

template <typename T>
class Array {
private:
    unsigned int count;
    unsigned int capacity;
    T* ptr;

public:
    Array() {
        count = 0;
        capacity = 0;
        ptr = nullptr;
    }

    Array(unsigned int initCount) {
        count = initCount;
        if (count == 0) {
            capacity = 0;
        }
        else if (count > 0x80000000) {
            throw "count exceeds max capacity";
        }
        else {
            capacity = 1;
            while (capacity < count) {
                capacity <<= 1;
            }
        }

        ptr = new T[capacity];
    }

    Array(unsigned int initCount, unsigned int initCapacity) {
        if (initCount > initCapacity) {
            throw "Initial count exceeds initial capacity";
        }

        count = initCount;
        capacity = initCapacity;
        ptr = new T[capacity];
    }

    Array(unsigned int initCount, T defaultValue) {
        count = initCount;
        if (count > 0x80000000) {
            throw "count exceeds max capacity";
        }
        else if (count == 0) {
            capacity = 0;
        }
        else {
            capacity = 1;
            while (capacity < count) {
                capacity <<= 1;
            }
        }

        ptr = new T[capacity]{ defaultValue };
    }

    Array(unsigned int initCount, unsigned int initCapacity, T defaultValue) {
        if (initCount > initCapacity) {
            throw "Initial count exceeds initial capacity";
        }

        count = initCount;
        capacity = initCapacity;
        ptr = new T[capacity]{ defaultValue };
    }

    Array(const Array& copy) {
        count = copy.count;
        capacity = copy.capacity;
        ptr = new T[capacity];
        for (unsigned int i = 0; i < count; i++) {
            ptr[i] = copy.ptr[i];
        }
    }

    ~Array() {
        delete[] ptr;
    }

    Array& operator=(const Array& copy) {
        delete[] ptr;

        count = copy.count;
        capacity = copy.capacity;
        ptr = new T[capacity];
        for (unsigned int i = 0; i < count; i++) {
            ptr[i] = copy.ptr[i];
        }

        return *this;
    }

    // viable if copy capacity is the name
    Array& quickCopy(const Array& copy) {
        count = copy.count;

        for (unsigned int i = 0; i < count; i++) {
            ptr[i] = copy.ptr[i];
        }

        return *this;
    }


    T& operator[](unsigned int index) {
        if (index < 0 || index > count - 1) {
            throw "array out of bounds";
        }

        return ptr[index];
    }

    const T& operator[](unsigned int index) const {
        if (index < 0 || index > count - 1) {
            throw "array out of bounds";
        }

        return ptr[index];
    }

    T* begin() {
        return ptr;
    }

    T* end() {
        return &ptr[count - 1] + 1;
    }

    unsigned int getCount() {
        return count;
    }

    unsigned int getCount() const {
        return count;
    }

    unsigned int getCapacity() {
        return capacity;
    }

    void setCapacity(unsigned int newCapacity) {
        T* oldPtr = ptr;

        capacity = newCapacity;
        ptr = new T[capacity];
        for (unsigned int i = 0; i < count; i++) {
            ptr[i] = oldPtr[i];
        }
        delete[] oldPtr;
    }

    // only works if type is valid to_string parameter
    std::string toString() {
        std::string output = "[";
        for (int i = 0; i < count; i++) {
            output += std::to_string(ptr[i]);
            if (i == count - 1) {
                output += "]";
                continue;
            }
            output += ",";
        }
        return output;
    }

    void resize() {
        T* oldPtr = ptr;

        if (capacity == 0) {
            capacity = 1;
        }

        capacity <<= 1;
        ptr = new T[capacity];
        for (unsigned int i = 0; i < count; i++) {
            ptr[i] = oldPtr[i];
        }
        delete[] oldPtr;
    }

    void append(T newObj) {
        if (count + 1 > capacity) {
            resize();
        }
        count++;
        ptr[count - 1] = newObj;
    }

    //void append(const Array<T>& arr) {
    //    int oldCount = count;
    //    count += arr.count;
    //    
    //    if (count > capacity) {
    //        capacity = 1;
    //        for (unsigned int i = 0; i < 32; i++) {
    //            if ((capacity << i) >= count) {
    //                capacity <<= i;
    //                break;
    //            }
    //        }

    //        T* oldPtr = ptr;
    //        ptr = new T[capacity];

    //        for (unsigned int i = 0; i < oldCount; i++) {
    //            ptr[i] = oldPtr[i];
    //        }
    //        delete[] oldPtr;

    //        for (unsigned int i = 0; i < arr.count; i++) {
    //            ptr[i + oldCount] = arr.ptr[i];
    //        }

    //        return;
    //    }

    //    for (unsigned int i = 0; i < arr.count; i++) {
    //        ptr[i + oldCount] = arr.ptr[i];
    //    }

    //    return;
    //}

    void remove(int index) {
        if (index < 0 || index > count - 1) {
            throw "array out of bounds";
        }

        count--;
        for (int i = index; i < count; i++) {
            ptr[i] = ptr[i + 1];
        }
        ptr[count] = T();
    }

    void resetCount() {
        count = 0;
    }

    void clear() {
        delete[] ptr;
        ptr = new T[capacity];
    }

    void clear(T defaultValue) {
        delete[] ptr;
        ptr = new T[capacity]{ defaultValue };
    }
};