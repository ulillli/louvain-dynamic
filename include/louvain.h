#pragma once
#include <map>
#include <set>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include "graph.h"

const double EPS = 1e-8;

template<class T>
class louvain {
	int N = 0;
	double M = 0;
	double first_modularity = 0.0;
	std::vector<int> teck_partition;
	std::vector<std::unordered_set<int>> teck_communities;
	std::vector<std::unordered_set<int>> communities;
	std::vector<int> result;

	graph<T> g;
	int teck_community_count = 0;
	std::vector<double> d;
	std::vector<std::unordered_map<int, double>> d_i_to_c; //номер сообщества, суммарный вес ребер-соседей, входящих в сообщество 
	std::vector<double> in;
	std::vector<double> tot;
	std::vector<double> in_G;
	std::vector<double> tot_G;
public:
	int getCommunitiesCount();
	std::vector<int> getPartition();
	std::vector<std::unordered_set<int>> getCommunities();
	double getFirstModularity();
	void printPartition();
	void printTeckPartition();
	void printInTot(int n);
	void printCommunities();
	void printDebugInfo();
	void printTeckCommunities();
	void printResultCommunities();
	void printD();

	void inizialization(const graph<T>& G);
	void reculculate(const graph<T>& g, const std::vector<int>& partition);
	void setSinglePartition(int n);
	double getModularity(const graph<T>& g, const std::vector<int>& partition);
	double getModularityOptimized(const graph<T>& g, const std::vector<int>& partition, bool flag);
	double getGain(const graph<T>& g, const int& v, const std::vector<int>& partition, const int& C);
	void remove(int v, int C, const graph<T>& g, std::vector<int>& partition);
	void insert(int v, int C, const graph<T>& g, std::vector<int>& partition);
	std::pair<double, int> getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition);
	void aggregateGraphOptimizedNew(graph<T>& g, std::vector<int>& partition);
	void aggregateGraphOptimizedNew2(graph<T>& g, std::vector<int>& partition);
	void aggregateGraphOptimized(graph<T>& g, std::vector<int>& partition);
	void moveNodes(graph<T>& g, std::vector<int>& partition);
	bool moveNodesNew(graph<T>& g, std::vector<int>& partition);
	louvain(const graph<T>& G);
};

template<class T>
int louvain<T>::getCommunitiesCount() { return teck_community_count; }
template<class T>
std::vector<std::unordered_set<int>> louvain<T>::getCommunities() { return communities; }
template<class T>
double louvain<T>::getFirstModularity() { return first_modularity; }
template<class T>
std::vector<int> louvain<T>::getPartition() {
	return result;
}
template<class T>
void louvain<T>::printPartition() {
	for (int i = 0; i < N; i++) std::cout << result[i] << " ";
	std::cout << std::endl;
}
template<class T>
void louvain<T>::printTeckPartition() {
	for (int i = 0; i < teck_community_count; i++) std::cout << teck_partition[i] << " ";
	std::cout << std::endl;
}
template<class T>
void louvain<T>::printInTot(int n) {
	std::cout << "in: ";
	for (int i = 0; i < n; i++) std::cout << in[i] << " ";
	std::cout << '\n';
	std::cout << "tot: ";
	for (int i = 0; i < n; i++) std::cout << tot[i] << " ";
	std::cout << '\n';
}
template<class T>
void louvain<T>::printCommunities() {
	for (int c = 0; c < teck_community_count; c++) {
		if (!communities[c].empty()) {
			std::cout << c << ": ";
			for (auto it = communities[c].begin(); it != communities[c].end(); it++) {
				int v = (*it);
				std::cout << v << " ";
			}
			std::cout << '\n';
		}
	}
}
template<class T>
void louvain<T>::printTeckCommunities() {
	
	for (int c = 0; c < teck_communities.size(); c++) {
		if (!teck_communities[c].empty()) {
			std::cout << "{";
			auto end = teck_communities[c].end();
			end--;
			for (auto it = teck_communities[c].begin(); it != end; it++) {
				int v = (*it);
				std::cout << v << ", ";
			}
			std::cout << *end << "}\n";
		}
	}
}
template<class T>
void louvain<T>::printResultCommunities() {
	std::cout << "[";
	for (int c = 0; c < teck_communities.size(); c++) {

		std::cout << "'";
		auto end = communities[c].end();
		end--;
		for (auto it = communities[c].begin(); it != end; it++) {
			int v = (*it);
			std::cout << v << ", ";
		}
		if (c == communities.size() - 1) std::cout << *end << "'";
		else std::cout << *end << "',";
	}
	std::cout << "]\n";
}
template<class T>
void louvain<T>::printDebugInfo() {
	std::cout << "partition: ";
	printTeckPartition();
	std::cout << "result: ";
	printPartition();
	printCommunities();
}
template<class T>
void louvain<T>::printD() {
	for (int i = 0; i < g.getVertexCount(); i++) std::cout << d[i] << " ";
	std::cout << std::endl;
}
template<class T>
void louvain<T>::reculculate(const graph<T>& g, const std::vector<int>& partition) {
	int n = g.getVertexCount();
	d = std::vector<double>(n, 0);
	d_i_to_c = std::vector<std::unordered_map<int, double>>(n);
	in = std::vector<double>(teck_community_count, 0);
	tot = std::vector<double>(teck_community_count, 0);
	for (int v = 0; v < n; v++) {
		int v_community = partition[v];
		in[v] = g.getWeightOfLoop(v);
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int neighbour_community = partition[neighbour];
			tot[v_community] += (*it).second;
			d[v] += (*it).second;
			if (v != neighbour) d_i_to_c[v][neighbour_community] += (*it).second;
		}
	}
}
template<class T>
void louvain<T>::setSinglePartition(int n) {
	teck_partition = std::vector<int>(teck_community_count);
	teck_communities = std::vector<std::unordered_set<int>>(teck_community_count);
	for (int i = 0; i < teck_community_count; i++) {
		teck_partition[i] = i;
		teck_communities[i].insert(i);
	}
	reculculate(g, teck_partition);
}
template<class T>
double louvain<T>::getGain(const graph<T>& g, const int& v, const std::vector<int>& partition, const int& C) {
	double totalWeight = 1.0 / (2.0 * g.getEdgeCount());
	double gain11 = (in[C] + 2.0 * d_i_to_c[v][C]) * totalWeight;
	double gain12 = ((tot[C] + d[v]) * totalWeight) * ((tot[C] + d[v]) * totalWeight);
	double gain21 = in[C] * totalWeight;
	double gain22 = (tot[C] * totalWeight) * (tot[C] * totalWeight);
	double gain23 = (d[v] * totalWeight) * (d[v] * totalWeight);
	return (gain11 - gain12) - (gain21 - gain22 - gain23);
}
template<class T>
void louvain<T>::remove(int v, int C, const graph<T>& g, std::vector<int>& partition) {
	in[C] = in[C] - 2.0 * d_i_to_c[v][C] - g.getWeightOfLoop(v);
	tot[C] = tot[C] - d[v];
	if (in[C] < 0) std::cout << "ERROR\n";
}
template<class T>
void louvain<T>::insert(int v, int C, const graph<T>& g, std::vector<int>& partition) {
	in[C] = in[C] + 2.0 * d_i_to_c[v][C] + g.getWeightOfLoop(v);
	tot[C] = tot[C] + d[v];
}
template<class T>
std::pair<double, int> louvain<T>::getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition) {
	int m = g.getEdgeCount();
	int v_community = partition[v];
	double tmp = in[v_community];
	remove(v, v_community, g, partition);
	double best_gain = getGain(g, v, partition, v_community);
	int best_community = v_community;
	for (auto it = d_i_to_c[v].begin(); it != d_i_to_c[v].end(); it++) {
		int community = it->first;
		if (community != v_community) {
			double gain = getGain(g, v, partition, community);
			if (gain > best_gain && gain > 0 ) { //добавлен gain > 0
				best_gain = gain;
				best_community = community;
			}
		}
	}
	insert(v, best_community, g, partition);
	return { best_gain,best_community };
}
template<class T>
void louvain<T>::aggregateGraphOptimized(graph<T>& g, std::vector<int>& partition) {
	//std::cout << "I'm in agregation\n";
	std::set<std::pair<int, int>> new_edges;
	int edges_count = 0;
	int n = g.getVertexCount();
	std::vector<std::unordered_set<int>> communities_ = std::vector<std::unordered_set<int>>();
	for (int c = 0; c < n; c++) {
		if (!teck_communities[c].empty()) communities_.push_back(teck_communities[c]);
	}
	teck_communities = communities_;
	if (teck_communities.size() == teck_community_count) return;
	communities_ = std::vector<std::unordered_set<int>>(teck_communities.size());
	for (int c_new = 0; c_new < teck_communities.size(); c_new++) {
		for (auto it = teck_communities[c_new].begin(); it != teck_communities[c_new].end(); it++) {
			int c_old = (*it); //это получается номер старого коммьюнити 
			partition[c_old] = c_new;
			for (auto it_ = communities[c_old].begin(); it_ != communities[c_old].end(); it_++) {
				int v = (*it_);
				result[v] = c_new;
				communities_[c_new].insert(v);
			}
		}
	}
	communities = communities_;
	teck_community_count = teck_communities.size();
	std::vector<std::vector<std::pair<int, double>>> adj(teck_community_count);
	std::vector<double> loops(teck_community_count, 0.0);
	for (int v = 0; v < n; v++) {

		int v_community = partition[v];
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int neighbour_community = partition[neighbour];
			bool flag = 0;
			for (int j = 0; j < adj[v_community].size(); j++) {
				if (adj[v_community][j].first == neighbour_community) {
					flag = 1;
					adj[v_community][j].second += (*it).second;
					break;
				}
			}
			if (!flag)adj[v_community].push_back({ neighbour_community,(*it).second });
			if (v_community == neighbour_community)  loops[v_community] += (*it).second;
			
		}
	}
	n = teck_community_count;
	double m = 0;
	for (int v = 0; v < n; v++) {
		for (auto it = adj[v].begin(); it != adj[v].end(); it++) {
			int u = it->first;
			double weight = it->second;
			m += weight;
		}
	}
	g = graph<T>(adj, loops, n, m / 2);
	setSinglePartition(teck_community_count);
}
template<class T>
void louvain<T>::aggregateGraphOptimizedNew(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();
	int count = 0;
	std::vector<std::unordered_set<int>> communities_;
	for (int c_old = 0; c_old < teck_communities.size(); c_old++) {
		if (!teck_communities[c_old].empty()) {
			communities_.push_back(std::unordered_set<int>());
			for (auto it = teck_communities[c_old].begin(); it != teck_communities[c_old].end(); it++) {
				int c = (*it); //это получается номер старого коммьюнити 
				partition[c] = count;
				for (auto it_ = communities[c].begin(); it_ != communities[c].end(); it_++) {
					int v = (*it_);
					result[v] = count;
					communities_[count].insert(v);
				}
			}
			count++;
		}
	}
	communities = communities_;
	teck_community_count = count;
	std::vector<std::vector<std::pair<int, double>>> adj(teck_community_count);
	std::vector<double> loops(teck_community_count, 0.0);
	std::vector<std::unordered_map<int, double>> edge_maps(teck_community_count);
	double m = 0;

	for (int v = 0; v < n; ++v) {
		int v_comm = partition[v];
		for (const auto& edge : g[v]) {
			int u = edge.first;
			double weight = edge.second;
			int u_comm = partition[u];
			m += weight;
			if (v_comm == u_comm)loops[v_comm] += weight;
			edge_maps[v_comm][u_comm] += weight;
		}
	}
	for (int i = 0; i < teck_community_count; ++i) {
		adj[i].reserve(edge_maps[i].size());
		for (const auto& [neighbor, weight] : edge_maps[i]) {
			adj[i].emplace_back(neighbor, weight);
		}
	}
	/*for (int v = 0; v < n; v++) {
		int v_community = partition[v];
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int neighbour_community = partition[neighbour];
			bool flag = 0;
			m += (*it).second;
			for (int j = 0; j < adj[v_community].size(); j++) {
				if (adj[v_community][j].first == neighbour_community) {
					flag = 1;
					
					adj[v_community][j].second += (*it).second;
					break;
				}
			}
			if (!flag) {
				
				adj[v_community].push_back({ neighbour_community,(*it).second });
			}
			if (v_community == neighbour_community)  loops[v_community] += (*it).second;

		}
	}*/
	n = teck_community_count;
	g = graph<T>(adj, loops, n, m / 2);
	setSinglePartition(teck_community_count);
}
template<class T>
void louvain<T>::aggregateGraphOptimizedNew2(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();

	std::vector<int> renumber(teck_community_count, -1);
	int new_count = 0;
	for (int v = 0; v < n; v++) {
		int c = partition[v];
		if (renumber[c] == -1) {
			renumber[c] = new_count++;
		}
		partition[v] = renumber[c];
	}

	if (new_count == n) return;

	for (int i = 0; i < N; i++) {
		result[i] = partition[result[i]];
	}
	teck_community_count = new_count;
	std::vector<std::vector<std::pair<int, double>>> adj(teck_community_count);
	std::vector<double> loops(teck_community_count, 0.0);
	std::vector<std::unordered_map<int, double>> edge_maps(teck_community_count);
	double m = 0;
	for (int v = 0; v < n; ++v) {
		int v_comm = partition[v];
		for (const auto& edge : g[v]) {
			int u = edge.first;
			double weight = edge.second;
			int u_comm = partition[u];
			m += weight;
			if (v_comm == u_comm)loops[v_comm] += weight;
			edge_maps[v_comm][u_comm] += weight;
		}
	}
	for (int i = 0; i < teck_community_count; ++i) {
		adj[i].reserve(edge_maps[i].size());
		for (const auto& [neighbor, weight] : edge_maps[i]) {
			adj[i].emplace_back(neighbor, weight);
		}
	}
	n = teck_community_count;
	g = graph<T>(adj, loops, n, m / 2);
	setSinglePartition(teck_community_count);
}
// в d_i_to_c учитываются связи только с другими вершинами, петли нет
template<class T>
void louvain<T>::moveNodes(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();
	double current_modularity = getModularityOptimized(g, partition, 1);
	double old_modularity = -1.0;
	std::cout << std::setprecision(15);
	do {
		old_modularity = current_modularity;
		for (int v = 0; v < n; v++) {
			int v_community = partition[v];
			std::pair<double, int> best_delta = getBestDelta(g, v, partition);
			if (best_delta.first > 0.0 && best_delta.second != v_community) {
				partition[v] = best_delta.second;
				teck_communities[best_delta.second].insert(v);
				teck_communities[v_community].erase(v);
				for (int j = 0; j < g[v].size(); j++) {
					if (v != g[v][j].first) {
						d_i_to_c[g[v][j].first][best_delta.second] += g[v][j].second;
						d_i_to_c[g[v][j].first][v_community] -= g[v][j].second;
						if (d_i_to_c[g[v][j].first][v_community] == 0) d_i_to_c[g[v][j].first].erase(v_community);
					}
				}
			}
		}
		//printTeckCommunities();
		current_modularity = getModularityOptimized(g, partition, 1);
	//std::cout << current_modularity << std::endl;
	} while (current_modularity > old_modularity);
}
template<class T>
bool louvain<T>::moveNodesNew(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();
	std::cout << std::setprecision(15);
	bool flag = 1;
	bool changed = 0;
	do {
		flag = 0;
		for (int v = 0; v < n; v++) {
			int v_community = partition[v];
			std::pair<double, int> best_delta = getBestDelta(g, v, partition);
			if (best_delta.first > 0.0 && best_delta.second != v_community) {
				flag = 1;
				changed = 1;
				partition[v] = best_delta.second;
				for (int j = 0; j < g[v].size(); j++) {
					if (v != g[v][j].first) {
						d_i_to_c[g[v][j].first][best_delta.second] += g[v][j].second;
						d_i_to_c[g[v][j].first][v_community] -= g[v][j].second;
						if (d_i_to_c[g[v][j].first][v_community] == 0) d_i_to_c[g[v][j].first].erase(v_community);
					}
				}
			}
		}
	} while (flag);
	return changed;
}

template<class T>
void louvain<T>::inizialization(const graph<T>& G) {
	N = G.getVertexCount();
	M = G.getEdgeCount();
	g = G;
	teck_community_count = N;
	teck_partition = std::vector<int>(N);
	result = std::vector<int>(N);

	setSinglePartition(teck_community_count);
	reculculate(g, teck_partition);
	communities = std::vector<std::unordered_set<int>>(teck_community_count);
	for (int i = 0; i < N; i++) {
		communities[i].insert(i);
		result[i] = i;
	}
	in_G = in;
	tot_G = tot;
}
template<class T>
louvain<T>::louvain(const graph<T>& G) {
	inizialization(G);
	bool flag = false;
	auto start = std::chrono::steady_clock::now();
	double old = getModularityOptimized(g, teck_partition, 0);
	first_modularity = old;
	std::cout << "Modularity: " << old << '\n';
	do {
		bool changed = moveNodesNew(g, teck_partition);
		//double curr = getModularityOptimized(g, teck_partition, 1);
		if (changed) {
			flag = true;
			aggregateGraphOptimizedNew2(g, teck_partition);
		}
		else flag = false;
	} while (flag);

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = end - start;
	std::cout << "Modularity: " << getModularity(G, result) << '\n';
	std::cout << "Count of communities: " << getCommunitiesCount() << std::endl;
	std::cout << time.count() << std::endl;
}
template<class T>
double louvain<T>::getModularityOptimized(const graph<T>& g, const std::vector<int>& partition, bool flag) {
	int community_count = 0;
	int n = g.getVertexCount(), m = g.getEdgeCount();
	community_count = teck_community_count;
	double result = 0;
	if (flag == 0) {
		for (int i = 0; i < community_count; i++) {
			result += in_G[i] / (2.0f * m) - (tot_G[i] / (2.0f * m)) * (tot_G[i] / (2.0f * m));
		}
	}
	else {
		for (int i = 0; i < community_count; i++) {
			result += in[i] / (2.0f * m) - (tot[i] / (2.0f * m)) * (tot[i] / (2.0f * m));
		}
	}
	return result;
}
template<class T>
double louvain<T>::getModularity(const graph<T>& g, const std::vector<int>& partition) {
	int community_count = 0;
	int n = g.getVertexCount(), m = g.getEdgeCount();
	for (int i = 0; i < n; i++) {
		if (partition[i] + 1 > community_count) community_count = partition[i] + 1;
	}
	std::vector<int> in = std::vector<int>(community_count, 0);
	std::vector<int> tot = std::vector<int>(community_count, 0);
	for (int v = 0; v < n; v++) {
		int v_community = partition[v];
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int neighbour_community = partition[neighbour];
			tot[v_community] += (*it).second;
			if (neighbour_community == v_community) {
				in[neighbour_community] += (*it).second;
			}
		}
	}

	double result = 0;
	for (int i = 0; i < community_count; i++) {
		result += in[i] / (2.0f * m) - (tot[i] / (2.0f * m)) * (tot[i] / (2.0f * m));
	}
	return result;
}
