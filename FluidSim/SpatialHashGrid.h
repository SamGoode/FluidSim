#pragma once
#include "raymath.h"
#include "Array.h"
#include "int2.h"

class SpatialHashGrid {
private:
    Vector2 m_pos;
    Vector2 m_size;
    int gridWidth;
    int gridHeight;
    float cellWidth;
    float cellHeight;
    int cellCount;
    int2 posOffsets[9];
    int hashOffsets[9];
    // (particle ID, cell ID/hash)
    Array<int2> hashList;
    Array<int> startIndices;
    Array<int> endIndices;
    Array<bool> dirty;
    Array<int> tempIDs;

public:
    SpatialHashGrid() {}
    SpatialHashGrid(Vector2 _pos, Vector2 _size, int _gridWidth, int _gridHeight);

    const Array<int2>& getHashList() { return hashList; }
    const Array<int>& getStartIndices() { return startIndices; }
    const Array<int>& getEndIndices() { return endIndices; }

    void setCellDirty(int cellHash) { dirty[cellHash] = true; }
    bool isCellDirty(int cellHash) { return dirty[cellHash]; }

    bool isValidCellPos(int2 cellPos) { return (cellPos.x >= 0 && cellPos.x < gridWidth && cellPos.y >= 0 && cellPos.y < gridHeight); }
    bool isValidPos(Vector2 pos) { return (pos.x >= m_pos.x && pos.x < m_pos.x + m_size.x && pos.y >= m_pos.y && pos.y < m_pos.y + m_size.y); }
    int2 getCellPos(int cellHash) { return { cellHash % gridWidth, cellHash / gridWidth }; }
    int2 getCellPos(Vector2 pos);

    int getCellHash(int2 cellPos) { return (cellPos.x + cellPos.y * gridWidth); }

    void generateHashList(const Array<Vector2> positions);
    void sortByCellHash();
    void generateLookup();

    const Array<int>& findWithin(int cellHash);
    const Array<int>& findNearby(int centreCellHash);
    const Array<int>& findNearby(int2 cellPos);

    void draw();
};