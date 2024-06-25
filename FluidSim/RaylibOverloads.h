#pragma once
#include "raymath.h"

static Vector2 operator+(Vector2 a, Vector2 b) {
    return { a.x + b.x, a.y + b.y };
}

static Vector2 operator-(Vector2 a, Vector2 b) {
    return { a.x - b.x, a.y - b.y };
}

static Vector2 operator*(Vector2 a, float scalar) {
    return { a.x * scalar, a.y * scalar };
}

static Vector2 operator/(Vector2 a, float scalar) {
    if (scalar == 0) {
        throw "divide by zero error";
    }
    return { a.x / scalar, a.y / scalar };
}

static void operator+=(Vector2& a, Vector2 b) {
    a = a + b;
}

static void operator-=(Vector2& a, Vector2 b) {
    a = a - b;
}

static bool operator==(Vector2 a, Vector2 b) {
    return (a.x == b.x && a.y == b.y);
}