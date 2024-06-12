#include "SpatialHashGrid.h"
#include <string>
#include "RadixSort.h"

SpatialHashGrid::SpatialHashGrid(Vector2 _pos, Vector2 _size, int _gridWidth, int _gridHeight) {
    m_pos = _pos;
    m_size = _size;
    gridWidth = _gridWidth;
    gridHeight = _gridHeight;
    cellWidth = (float)(m_size.x / gridWidth);
    cellHeight = (float)(m_size.y / gridHeight);
    cellCount = gridWidth * gridHeight;

    for (int i = 0; i < 9; i++) {
        posOffsets[i] = { (i % 3) - 1, (i / 3) - 1 };
        hashOffsets[i] = { ((i % 3) - 1) + ((i / 3) - 1) * gridWidth };
    }

    hashList = Array<int2>(0);
    startIndices = Array<int>(cellCount);
    endIndices = Array<int>(cellCount);
    dirty = Array<bool>(cellCount);
    tempIDs = Array<int>(0);
}

int2 SpatialHashGrid::getCellPos(Vector2 pos) {
    Vector2 adjusted = Vector2Subtract(pos, m_pos);
    int2 cellPos = { (int)floor(adjusted.x / cellWidth), (int)floor(adjusted.y / cellHeight) };
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

void SpatialHashGrid::generateHashList(const Array<Vector2> positions) {
    hashList = Array<int2>(positions.getCount());

    for (int i = 0; i < hashList.getCount(); i++) {
        hashList[i] = { i, getCellHash(getCellPos(positions[i]))};
    }

    tempIDs = Array<int>(0, hashList.getCapacity());
}

void SpatialHashGrid::sortByCellHash() {
    if (hashList.getCount() <= 1) {
        return;
    }

    //std::sort(&hashList[0], &hashList[hashList.getCount() - 1] + 1, compareCellHash);
    radixSort(hashList, 127);
}

void SpatialHashGrid::generateLookup() {
    startIndices.clear(-1);
    endIndices.clear(-1);
    dirty.clear(false);

    int currentStart = 0;
    int previousCellHash = -1;
    for (int i = 0; i < hashList.getCount() - 1; i++) {
        if (hashList[i].y == hashList[i + 1].y) {
            continue;
        }

        startIndices[hashList[i].y] = currentStart;

        if (previousCellHash != -1) {
            endIndices[previousCellHash] = currentStart - 1;
        }
        previousCellHash = hashList[i].y;

        currentStart = i + 1;
    }

    endIndices[previousCellHash] = currentStart - 1;
    startIndices[hashList[hashList.getCount() - 1].y] = currentStart;
    endIndices[hashList[hashList.getCount() - 1].y] = hashList.getCount() - 1;
}

const Array<int>& SpatialHashGrid::findWithin(int cellHash) {
    tempIDs.resetCount();
    int startIndex = startIndices[cellHash];
    int endIndex = endIndices[cellHash];

    if (startIndex < 0 || endIndex < 0) {
        return tempIDs;
    }

    for (int i = startIndex; i < endIndex + 1; i++) {
        tempIDs.append(hashList[i].x);
    }

    return tempIDs;
}

const Array<int>& SpatialHashGrid::findNearby(int centreCellHash) {
    tempIDs.resetCount();
    int2 centreCellPos = getCellPos(centreCellHash);
    int cellHash;
    int startIndex;
    int endIndex;

    for (int i = 0; i < 9; i++) {
        if (!isValidCellPos(centreCellPos + posOffsets[i])) {
            continue;
        }

        cellHash = centreCellHash + hashOffsets[i];

        if (dirty[cellHash]) {
            continue;
        }

        startIndex = startIndices[cellHash];
        endIndex = endIndices[cellHash];

        if (startIndex < 0 || endIndex < 0) {
            continue;
        }

        for (int n = startIndex; n < endIndex + 1; n++) {
            tempIDs.append(hashList[n].x);
        }
    }

    return tempIDs;
}

const Array<int>& SpatialHashGrid::findNearby(int2 cellPos) {
    tempIDs.resetCount();
    int centreCellHash = getCellHash(cellPos);

    for (int i = 0; i < 9; i++) {
        if (!isValidCellPos(cellPos + posOffsets[i])) {
            continue;
        }

        int cellHash = centreCellHash + hashOffsets[i];

        int startIndex = startIndices[cellHash];
        int endIndex = endIndices[cellHash];
        if (startIndex < 0 || endIndex < 0) {
            continue;
        }

        for (int n = startIndex; n < endIndex + 1; n++) {
            tempIDs.append(hashList[n].x);
        }
    }

    return tempIDs;
}

void SpatialHashGrid::draw() {
    Vector2 pos = GetMousePosition();
    

    int2 cellPos = getCellPos(pos);
    Array<int> withinIDs = findWithin(getCellHash(cellPos));
    Array<int> IDs = findNearby(cellPos);

    int x = cellPos.x;
    int y = cellPos.y;
    DrawRectangle(m_pos.x + (x - 1) * cellWidth, m_pos.y + (y - 1) * cellHeight, cellWidth * 3, cellHeight * 3, ORANGE);
    DrawRectangle(m_pos.x + x * cellWidth, m_pos.y + y * cellHeight, cellWidth, cellHeight, RED);

    for (int i = 0; i < gridWidth; i++) {
        for (int n = 0; n < gridHeight; n++) {
            DrawRectangleLines(round(m_pos.x + i * cellWidth), round(m_pos.y + n * cellHeight), cellWidth, cellHeight, BLACK);
        }
    }

    DrawRectangleLines(m_pos.x, m_pos.y, m_size.x, m_size.y, BLUE);

    //for (int i = 0; i < hashList.getCount(); i++) {
    //    bool isNearby = false;
    //    for (int n = 0; n < IDs.getCount(); n++) {
    //        if (hashList[i].x == IDs[n]) {
    //            isNearby = true;
    //        }
    //    }

    //    Vector2 c_pos = critters[hashList[i].x].GetPosition();

    //    if (isNearby) {
    //        DrawCircle(c_pos.x, c_pos.y, 15, RED);
    //    }
    //    else {
    //        DrawCircle(c_pos.x, c_pos.y, 15, GREEN);
    //    }
    //}

    //DrawCircle(pos.x, pos.y, 15, BLUE);

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