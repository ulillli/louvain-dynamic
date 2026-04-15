#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "graph.h"
#define LOOP_TIME 0
struct Edge {
    int from;       
    int to;        
    double weight;    
    long long time;
};

class DynamicGraph {
private:
    std::vector<std::vector<Edge>> adj;
    std::vector<double> loops;

    int max_vertex_id;
    long long min_time;
    long long max_time;

public:
    DynamicGraph() : max_vertex_id(0), min_time(std::numeric_limits<long long>::max()), max_time(0) {}
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "File " << filename << " isn't opened" << std::endl;
            return false;
        }

        std::string line;
        std::vector<Edge> temp_edges;
        max_vertex_id = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            int u, v;
            long long t;
            double w = 1; 
            if (iss >> u >> v >> t && u!=v) {
                temp_edges.push_back({ u, v, w, t });
                max_vertex_id = std::max({ max_vertex_id, u, v });
                min_time = std::min(min_time, t);
                max_time = std::max(max_time, t);
            }
        }
        file.close();
        adj.assign(max_vertex_id + 1, std::vector<Edge>());
        loops.assign(max_vertex_id + 1,0.0);
        for (const auto& edge : temp_edges) {
            adj[edge.from].push_back(edge);
            adj[edge.to].push_back({edge.to, edge.from,edge.weight,edge.time});
        }
        for (int i = 0; i < max_vertex_id+1; i++) { //добавляем петли сами с 0-ым весом
            adj[i].push_back({i,i,loops[i],min_time});
        }
        return true;
    }
    long long getMinTime() const { return min_time; }
    long long getMaxTime() const { return max_time; }
    int getMaxVertexId() const { return max_vertex_id; }
    graph<std::vector<std::pair<int, double>>> getSnapshot(long long t_start, long long t_end) const {
        int N = max_vertex_id + 1;
        std::vector<std::vector<std::pair<int, double>>> snap_adj(N);
        double M_count = 0;
        for (int u = 0; u < N; ++u) {
            //std::unordered_map<int, double> aggregated_edges; // используем map для суммирования весов, если в этом окне вершина u обратилась к v несколько раз.
            //for (const auto& edge : adj[u]) {
            //    if (edge.time >= t_start && edge.time <= t_end) {
            //        if (edge.from == edge.to) {
            //            snap_loops[edge.from] = edge.weight;
            //            aggregated_edges[edge.to] = edge.weight;
            //        }
            //        else {
            //            aggregated_edges[edge.to] = edge.weight;
            //            aggregated_edges[edge.from] = edge.weight;
            //        }
            //    }
            //}
            //for (const auto& pair : aggregated_edges) {
            //    snap_adj[u].push_back({ pair.first, pair.second });
            //    M_count += pair.second;
            //}
            //M_count += snap_loops[u]; 
            for (const auto& edge : adj[u]) {
                //if (edge.time>=t_start && edge.time <= t_end) {
                if (edge.time <= t_end) {
                    M_count += edge.weight;
                    snap_adj[edge.from].push_back({ edge.to,edge.weight });
                }
            }
        }
        
        return  graph<std::vector<std::pair<int, double>>>(snap_adj, loops, N, M_count/2);
    }

    void printGraph() const {
        std::cout << "--- Main graph (adj) ---" << std::endl;
        for (int i = 0; i <= max_vertex_id; ++i) {
            if (!adj[i].empty()) {
                std::cout << "Вершина " << i << ":\n";
                for (const auto& e : adj[i]) {
                    std::cout << "  -> " << e.to << " (Вес: " << e.weight
                        << ", Время: " << e.time << ")\n";
                }
            }
        }
        std::cout << "\n--- Только петли (loops) ---" << std::endl;
        for (int i = 0; i <= max_vertex_id; ++i) {
            std::cout << loops[i] << " ";
        }
        std::cout << std::endl;
    }
};
