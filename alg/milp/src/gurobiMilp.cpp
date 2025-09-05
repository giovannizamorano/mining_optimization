#include <bits/stdc++.h>
#include "gurobi_c++.h"

struct Block {
  int64_t id;
  double X, Y, Z;
  double density;
  double economic; // valor por bloque
  double value;    // = economic
  int ix, iy, iz;
};

struct Key {
  int ix, iy, iz;
  bool operator==(const Key& o) const { return ix==o.ix && iy==o.iy && iz==o.iz; }
};
struct KeyHash {
  std::size_t operator()(const Key& k) const noexcept {
    return ( (std::hash<int>()(k.ix) * 73856093u)
           ^ (std::hash<int>()(k.iy) * 19349663u)
           ^ (std::hash<int>()(k.iz) * 83492791u) );
  }
};

static std::vector<std::string> splitCSV(const std::string& line) {
  std::vector<std::string> out; std::string cur; bool inQ=false;
  for (char c: line) {
    if (c=='"') inQ=!inQ;
    else if (c==',' && !inQ) { out.push_back(cur); cur.clear(); }
    else cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "Uso: " << argv[0] << " <ruta/MARVIN_BM.csv> [--limit=N]\n";
      return 1;
    }
    std::string csvPath = argv[1];
    int64_t limitN = -1;
    for (int i=2;i<argc;++i) {
      std::string a = argv[i];
      if (a.rfind("--limit=",0)==0) limitN = std::stoll(a.substr(8));
      else { std::cerr << "Arg no reconocido: " << a << "\n"; return 1; }
    }

    // Abrir CSV
    std::ifstream fin(csvPath);
    if (!fin) { std::cerr << "No se pudo abrir " << csvPath << "\n"; return 1; }
    std::string header; std::getline(fin, header);
    auto cols = splitCSV(header);
    auto findIdx = [&](const std::string& name)->int {
      for (int i=0;i<(int)cols.size();++i) if (cols[i]==name) return i;
      return -1;
    };
    int ixX=findIdx("X"), ixY=findIdx("Y"), ixZ=findIdx("Z");
    int ixDensity=findIdx("Density"), ixEconomic=findIdx("Economic");
    if (ixX<0||ixY<0||ixZ<0||ixDensity<0||ixEconomic<0) {
      std::cerr << "Faltan encabezados. Se requieren: X,Y,Z,Density,Economic\n";
      return 1;
    }

    std::vector<Block> blocks; blocks.reserve(70000);
    std::vector<double> allX, allY, allZ; allX.reserve(1000); allY.reserve(1000); allZ.reserve(1000);
    std::string line; int64_t id=0;
    while (std::getline(fin, line)) {
      if (line.empty()) continue;
      auto v = splitCSV(line);
      if ((int)v.size() <= std::max({ixX,ixY,ixZ,ixDensity,ixEconomic})) continue;
      Block b;
      b.id = id++;
      b.X  = std::stod(v[ixX]);
      b.Y  = std::stod(v[ixY]);
      b.Z  = std::stod(v[ixZ]);
      b.density  = std::stod(v[ixDensity]);
      b.economic = std::stod(v[ixEconomic]);
      b.value    = b.economic; // ya viene por bloque
      blocks.push_back(b);
      allX.push_back(b.X); allY.push_back(b.Y); allZ.push_back(b.Z);
      if (limitN>0 && (int64_t)blocks.size()>=limitN) break;
    }
    fin.close();
    if (blocks.empty()) { std::cerr << "CSV vacío.\n"; return 1; }

    auto uniq = [](std::vector<double>& a){
      std::sort(a.begin(), a.end()); a.erase(std::unique(a.begin(), a.end()), a.end());
    };
    uniq(allX); uniq(allY); uniq(allZ);
    if (allX.size()<2||allY.size()<2||allZ.size()<2) {
      std::cerr << "Insuficiente muestreo para inferir malla.\n"; return 1;
    }
    double dx = allX[1]-allX[0], dy = allY[1]-allY[0], dz = allZ[1]-allZ[0];
    std::cerr << "Malla inferida: dx="<<dx<<" dy="<<dy<<" dz="<<dz<<" (m)\n";

    // Mapear a índices de rejilla
    double minX=allX.front(), minY=allY.front(), minZ=allZ.front();
    auto toIndex = [&](double c, double m, double d)->int {
      return (int)std::llround((c-m)/d);
    };
    std::unordered_map<Key,int64_t,KeyHash> idAt; idAt.reserve(blocks.size()*1.3);
    for (auto& b: blocks) {
      b.ix = toIndex(b.X,minX,dx); b.iy=toIndex(b.Y,minY,dy); b.iz=toIndex(b.Z,minZ,dz);
      idAt[{b.ix,b.iy,b.iz}] = b.id;
    }

    // Talud 45° => R=1 (3x3 arriba)
    int R = 1;

    // Construir arcos de precedencia (i -> j superior)
    std::vector<std::pair<int64_t,int64_t>> arcs; arcs.reserve(blocks.size()*9);
    for (const auto& b: blocks) {
      Key k; k.iz = b.iz + 1;
      for (int u=-R; u<=R; ++u)
        for (int v=-R; v<=R; ++v) {
          k.ix = b.ix + u; k.iy = b.iy + v;
          auto it = idAt.find(k);
          if (it!=idAt.end()) arcs.emplace_back(b.id, it->second);
        }
    }
    std::cerr << "Arcos de precedencia: " << arcs.size() << "\n";

    // Modelo Gurobi
    GRBEnv env = GRBEnv(true);
    env.set("LogToConsole","1");
    env.start();
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "UPL_single_period_MILP");

    // Variables x_b (binarias). Density<=0 => x_b fijada a 0
    std::vector<GRBVar> x; x.reserve(blocks.size());
    for (const auto& b: blocks) {
      double lb=0.0, ub=(b.density<=0.0?0.0:1.0);
      x.push_back(model.addVar(lb, ub, b.value, GRB_BINARY, "x_"+std::to_string(b.id)));
    }

    // Objetivo
    GRBLinExpr obj = 0.0;
    for (size_t i=0;i<blocks.size();++i) obj += blocks[i].value * x[i];
    model.setObjective(obj, GRB_MAXIMIZE);

    // Precedencias: x_i <= x_j
    for (auto &e: arcs) model.addConstr(x[e.first] <= x[e.second], "prec");

    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if (status==GRB_OPTIMAL || status==GRB_SUBOPTIMAL) {
      double z = model.get(GRB_DoubleAttr_ObjVal);
      std::cout << std::fixed << std::setprecision(2);
      std::cout << "Objective: " << z << "\n";
      std::ofstream fout("pit_solution.csv");
      fout << "X,Y,Z,value,chosen\n";
      for (size_t i=0;i<blocks.size();++i) {
        int chosen = (x[i].get(GRB_DoubleAttr_X) > 0.5) ? 1 : 0;
        fout << blocks[i].X << "," << blocks[i].Y << "," << blocks[i].Z << ","
             << blocks[i].value << "," << chosen << "\n";
      }
      fout.close();
      std::cout << "Solución escrita en pit_solution.csv\n";
    } else {
      std::cerr << "Optimización no óptima. Status="<<status<<"\n";
    }
  } catch (GRBException &e) {
    std::cerr << "Gurobi error: " << e.getMessage() << " code=" << e.getErrorCode() << "\n";
    return 1;
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
