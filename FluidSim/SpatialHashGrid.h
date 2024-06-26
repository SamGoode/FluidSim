#pragma once
#include "raylib.h"
#include "Array.h"
#include "int2.h"
#include "Particle.h"

class SpatialHashGrid {
private:
    Vector2 size;
    int gridWidth;
    int gridHeight;
    float cellWidth;
    float cellHeight;
    int cellCount;

    Array<int> hashOffsets;

    // (particle ID, cell ID/hash)
    Array<int2> hashList;
    // (start Index, end Index)
    Array<int2> indexLookup;

public:
    SpatialHashGrid() {}
    SpatialHashGrid(Vector2 _size, int _gridWidth, int _gridHeight);
    SpatialHashGrid(Vector2 _size, float _cellWidth, float _cellHeight);

    const Array<int>& getHashOffsets() { return hashOffsets; }

    const Array<int2>& getHashList() { return hashList; }
    const Array<int2>& getIndexLookup() { return indexLookup; }

    bool isValidCellPos(int2 cellPos) { return (cellPos.x >= 0 && cellPos.x < gridWidth && cellPos.y >= 0 && cellPos.y < gridHeight); }
    int2 getCellPos(int cellHash) { return { cellHash % gridWidth, cellHash / gridWidth }; }
    int2 getCellPos(Vector2 pos);
    bool isValidCellHash(int cellHash) { return { cellHash >= 0 && cellHash < cellCount }; }
    int getCellHash(int2 cellPos) { return (cellPos.x + cellPos.y * gridWidth); }

    void generateHashList(const Array<Particle>& particles, const Array<int>& objectPool, int activeCount);
    void sortByCellHash();
    void generateLookup();

    void draw(Vector2 pos, float scale, Vector2 testPos);
};