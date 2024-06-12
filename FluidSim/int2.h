#pragma once

struct int2 {
public:
    int x;
    int y;

public:
    friend int2 operator+(int2 a, int2 b) {
        return { a.x + b.x, a.y + b.y };
    }
};