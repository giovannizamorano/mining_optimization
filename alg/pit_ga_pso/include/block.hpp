#pragma once
#include <vector>
#include <unordered_map>

struct Block {
    int     id;
    double  value;
    double  tonnage;
    double  cu_grade;    // nueva
    double  au_grade;    // nueva
};

using BlockArray = std::vector<Block>;
using ArcAdj     = std::unordered_map<int, std::vector<int>>; // pred[v] = {parents}
