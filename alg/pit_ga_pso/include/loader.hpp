#pragma once
#include <string>
#include "block.hpp"
#include <vector>
#include <unordered_map>
#include <nlohmann/json.hpp>

// Carga bloques y precedencias desde JSON generado por data_generator.py
BlockArray load_blocks_json(const std::string &json_path);
ArcAdj     load_arcs_json  (const std::string &json_path);
