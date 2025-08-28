#include "loader.hpp"
#include <fstream>
#include <stdexcept>

BlockArray load_blocks_json(const std::string &fn) {
    std::ifstream f(fn);
    if (!f) throw std::runtime_error("No se pudo abrir JSON de entrada: " + fn);
    nlohmann::json j; 
    f >> j;
    const auto &jblocks = j.at("blocks");
    BlockArray B;
    B.reserve(jblocks.size() + 1);
    B.push_back({0, 0.0, 0.0, 0.0, 0.0}); // dummy para indexación 1-based
    for (const auto &bj : jblocks) {
        Block b;
        b.id       = bj.at("id").get<int>();
        b.value    = bj.at("value").get<double>();
        b.tonnage  = bj.at("tonnage").get<double>();
        // Si tienes cu/au en tu Block struct, inclúyelos:
        b.cu_grade = bj.at("cu").get<double>();
        b.au_grade = bj.at("au").get<double>();
        B.push_back(b);
    }
    return B;
}

ArcAdj load_arcs_json(const std::string &fn) {
    std::ifstream f(fn);
    if (!f) throw std::runtime_error("No se pudo abrir JSON de entrada: " + fn);
    nlohmann::json j;
    f >> j;
    const auto &jprec = j.at("precedences");
    ArcAdj pred;
    for (auto it = jprec.begin(); it != jprec.end(); ++it) {
        int v = std::stoi(it.key());
        pred[v] = it.value().get<std::vector<int>>();
    }
    return pred;
}
