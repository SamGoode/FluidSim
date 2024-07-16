#include "SpatialHashGrid.h"
#include <string>
#include "RadixSort.h"

SpatialHashGrid::SpatialHashGrid(Vector2 _size, int _gridWidth, int _gridHeight) {
    size = _size;
    gridWidth = _gridWidth;
    gridHeight = _gridHeight;
    cellWidth = (float)(size.x / _gridWidth);
    cellHeight = (float)(size.y / _gridHeight);
    cellCount = _gridWidth * _gridHeight;

    // hash offset values for 3x3 grid
    hashOffsets = Array<int>(9);
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            hashOffsets[x + y * 3] = (x - 1) + (y - 1) * gridWidth;
        }
    }

    hashList = Array<int2>(0);
    indexLookup = Array<int2>(cellCount);
}

SpatialHashGrid::SpatialHashGrid(Vector2 _size, float _cellWidth, float _cellHeight) {
    size = _size;
    gridWidth = size.x / _cellWidth;
    gridHeight = size.y / _cellHeight;
    cellWidth = _cellWidth;
    cellHeight = _cellHeight;
    cellCount = gridWidth * gridHeight;

    // hash offset values for 3x3 grid
    hashOffsets = Array<int>(9);
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            hashOffsets[x + y * 3] = (x - 1) + (y - 1) * gridWidth;
        }
    }

    hashList = Array<int2>(0);
    indexLookup = Array<int2>(cellCount);
}

int2 SpatialHashGrid::getCellPos(Vector2 pos) {
    int2 cellPos = { (int)floor(pos.x / cellWidth), (int)floor(pos.y / cellHeight) };
    if (cellPos.x < 0) {
        cellPos.x = 0;
    }
    else if (cellPos.x >= gridWidth) {
        cellPos.x = gridWidth - 1;
    }

    if (cellPos.y < 0) {
        cellPos.y = 0;
    }
    else if (cellPos.y >= gridHeight) {
        cellPos.y = gridHeight - 1;
    }

    return cellPos;
}

void SpatialHashGrid::generateHashList(const Array<Vector2>& positions, const Array<int>& objectPool, int activeCount) {
    hashList = Array<int2>(activeCount);
    
    for (int i = 0; i < activeCount; i++) {
        int particleID = objectPool[i];
        hashList[i] = { particleID, getCellHash(getCellPos(positions[particleID])) };
    }
}

void SpatialHashGrid::sortByCellHash() {
    if (hashList.getCount() <= 1) {
        return;
    }

    radixSort(hashList, cellCount - 1);
}

void SpatialHashGrid::generateLookup() {
    indexLookup.clear({ -1, -1 });

    if (hashList.getCount() == 0) {
        return;
    }

    int currentStart = 0;
    int previousCellHash = -1;
    for (int i = 0; i < hashList.getCount() - 1; i++) {
        if (hashList[i].y == hashList[i + 1].y) {
            continue;
        }

        indexLookup[hashList[i].y].x = currentStart;

        if (previousCellHash != -1) {
            indexLookup[previousCellHash].y = currentStart - 1;
        }
        previousCellHash = hashList[i].y;

        currentStart = i + 1;
    }

    if (previousCellHash != -1) {
        indexLookup[previousCellHash].y = currentStart - 1;
    }

    indexLookup[hashList[hashList.getCount() - 1].y].x = currentStart;
    indexLookup[hashList[hashList.getCount() - 1].y].y = hashList.getCount() - 1;
}

void SpatialHashGrid::draw(Vector2 pos, float scale, Vector2 testPos) {
    int2 cellPos = getCellPos(testPos);

    int x = cellPos.x;
    int y = cellPos.y;
    DrawRectangle(pos.x +(x - 1) * cellWidth * scale, pos.y + (y - 1) * cellHeight * scale, cellWidth * scale * 3, cellHeight * scale * 3, ORANGE);
    DrawRectangle(pos.x + x * cellWidth * scale, pos.y + y * cellHeight * scale, cellWidth * scale, cellHeight * scale, RED);

    for (int i = 1; i < gridWidth; i++) {
        DrawLine(pos.x + i * cellWidth * scale, pos.y, pos.x + i * cellWidth * scale, pos.y + size.y * scale, BLACK);
    }
    for (int i = 1; i < gridHeight; i++) {
        DrawLine(pos.x, pos.y + i * cellHeight * scale, pos.x + size.x * scale, pos.y + i * cellHeight * scale, BLACK);
    }

    DrawRectangleLines(pos.x, pos.y, size.x * scale, size.y * scale, BLACK);

    //std::string cellID = "Cell ID: " + std::to_string(getCellHash(cellPos));
    //DrawText(cellID.c_str(), 10, 100, 20, PURPLE);

    //std::string nearby = "Within bounds: " + std::to_string(IDs.getCount());
    //DrawText(nearby.c_str(), 10, 140, 20, PURPLE);

    //std::string within = "Within cell: " + std::to_string(withinIDs.getCount());
    //DrawText(within.c_str(), 10, 160, 20, PURPLE);

    //std::string valid = isValidPos(GetMousePosition()) ? "true" : "false";
    //valid = "Valid location: " + valid;
    //DrawText(valid.c_str(), 10, 180, 20, PURPLE);

    //std::string IDsText = "IDs: [";
    //for (int i = 0; i < IDs.getCount(); i++) {
    //    IDsText += std::to_string(IDs[i]);
    //    if (i != IDs.getCount() - 1) {
    //        IDsText += ",";
    //    }
    //}
    //IDsText += "]";
    //DrawText(IDsText.c_str(), 200, 20, 20, PURPLE);
}