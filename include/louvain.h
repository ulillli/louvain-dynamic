#pragma once
#include <map>
#include <set>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include "graph.h"

#include <algorithm>
#include <random>
#include <numeric>
#include <omp.h>

const double EPS = 1e-8;

template<class T>
class louvain {
	int N = 0;
	double M = 0;
	double first_modularity = 0.0;
	std::vector<int> teck_partition;
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
	double getFirstModularity();
	void printPartition();
	void printTeckPartition();
	void printInTot(int n);
	void printDebugInfo();
	void printD();

	void inizialization(const graph<T>& G);
	void reculculate(const graph<T>& g, const std::vector<int>& partition);
	void setSinglePartition(int n);
	double getModularity(const graph<T>& g, const std::vector<int>& partition);
	double getModularityOptimized(const graph<T>& g, const std::vector<int>& partition, bool flag);
	double getGain(const graph<T>& g, const int& v, const int& C, double weight_to_c, double totalWeight);
	void remove(int v, int C, const graph<T>& g, double weight_to_c);
	void insert(int v, int C, const graph<T>& g, double weight_to_c);
	bool getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition);
	void aggregateGraph(graph<T>& g, std::vector<int>& partition);
	bool moveNodes(graph<T>& g, std::vector<int>& partition);
	bool moveNodesRandom(graph<T>& g, std::vector<int>& partition);
	louvain(const graph<T>& G);
};

template<class T>
int louvain<T>::getCommunitiesCount() { return teck_community_count; }
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
void louvain<T>::printDebugInfo() {
	std::cout << "partition: ";
	printTeckPartition();
	std::cout << "result: ";
	printPartition();
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
		double loop = g.getWeightOfLoop(v);
		int v_community = partition[v];
		in[v] = loop;
		tot[v_community] += loop;
		d[v] += loop;
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int n_community = partition[neighbour];
			double weight = (*it).second;

			tot[v_community] += weight;
			d[v] += weight;
			d_i_to_c[v][n_community] += weight;
		}
	}
}
template<class T>
void louvain<T>::setSinglePartition(int n) {
	teck_partition = std::vector<int>(teck_community_count);
	for (int i = 0; i < teck_community_count; i++) {
		teck_partition[i] = i;
	}
	reculculate(g, teck_partition);
}
template<class T>
double louvain<T>::getGain(const graph<T>& g, const int& v, const int& C, double weight_to_c, double totalWeight) {
	double gain11 = (in[C] + 2.0 * weight_to_c) * totalWeight;
	double gain12 = ((tot[C] + d[v]) * totalWeight) * ((tot[C] + d[v]) * totalWeight);
	double gain21 = in[C] * totalWeight;
	double gain22 = (tot[C] * totalWeight) * (tot[C] * totalWeight);
	double gain23 = (d[v] * totalWeight) * (d[v] * totalWeight);
	return (gain11 - gain12) - (gain21 - gain22 - gain23);
}
template<class T>
void louvain<T>::remove(int v, int C, const graph<T>& g, double weight_to_c) {
	in[C] = in[C] - 2.0 * weight_to_c - g.getWeightOfLoop(v);
	tot[C] = tot[C] - d[v];
	if (in[C] < 0.0) std::cout << "ERROR: " << in[C] ;
}
template<class T>
void louvain<T>::insert(int v, int C, const graph<T>& g, double weight_to_c) {
	in[C] = in[C] + 2.0 * weight_to_c + g.getWeightOfLoop(v);
	tot[C] = tot[C] + d[v];
}
template<class T>
bool louvain<T>::getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition) {
	bool flag = false;
	double totalWeight = 1.0 / (2.0 * g.getEdgeCount());

	int v_community = partition[v];
	double tmp = in[v_community];

	double weight_to_old = 0.0;
	auto it_old = d_i_to_c[v].find(v_community);
	if (it_old != d_i_to_c[v].end()) {
		weight_to_old = it_old->second;
	}

	remove(v, v_community, g, weight_to_old);
	double best_gain = getGain(g, v, v_community, weight_to_old, totalWeight);
	int best_community = v_community;

	for (auto it = d_i_to_c[v].begin(); it != d_i_to_c[v].end(); it++) {
		int community = it->first;
		if (community != v_community) {
			double gain = getGain(g, v, community, it->second, totalWeight);
			if (gain > best_gain && gain > 0 ) {
				best_gain = gain;
				best_community = community;
			}
		}
	}
	double weight_to_best = weight_to_old;
	if (best_community != v_community) {
		auto it_best = d_i_to_c[v].find(best_community);
		if (it_best != d_i_to_c[v].end()) {
			weight_to_best = it_best->second;
		}
	}
	insert(v, best_community, g, weight_to_best);
	if (best_gain > 0.0 && best_community != v_community) {
		flag = true;
		partition[v] = best_community;
		for (int j = 0; j < g[v].size(); j++) {
			int neighbor = g[v][j].first;
			double edge_w = g[v][j].second;
			d_i_to_c[neighbor][best_community] += edge_w;
			d_i_to_c[neighbor][v_community] -= edge_w;
			if (d_i_to_c[neighbor][v_community] <= 0.0) {
				d_i_to_c[neighbor].erase(v_community);
			}
			
		}
	}
	return flag;
}
template<class T>
void louvain<T>::aggregateGraph(graph<T>& g, std::vector<int>& partition) {
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
		loops[v_comm] += g.getWeightOfLoop(v);
		for (const auto& edge : g[v]) {
			int u = edge.first;
			double weight = edge.second;
			int u_comm = partition[u];
			m += weight;
			if (v_comm == u_comm)loops[v_comm] += weight;
			else edge_maps[v_comm][u_comm] += weight;
		}
	}
	for (int i = 0; i < teck_community_count; ++i) {
		adj[i].reserve(edge_maps[i].size());
		for (const auto& [neighbor, weight] : edge_maps[i]) {
			adj[i].emplace_back(neighbor, weight);
		}
	}
	n = teck_community_count;
	g = graph<T>(adj, loops, n, this->M);
	setSinglePartition(teck_community_count);
	std::cout << "New graph with n = " << n << " and communities count = " << teck_community_count << std::endl;
}
template<class T>
bool louvain<T>::moveNodes(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();

	bool global_changed = false;
	bool pass_changed = false;

	do {
		pass_changed = false;
		for (int v = 0; v < n; v++) {

			bool moved = getBestDelta(g, v, partition);

			if (moved) {
				pass_changed = true;
			}
		}
		if (pass_changed) {
			global_changed = true;
		}

	} while (pass_changed); 

	return global_changed;
}
template<class T>
bool louvain<T>::moveNodesRandom(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();
	std::cout << std::setprecision(15);

	bool changed_global = false; 
	bool changed_in_pass = false;

	std::vector<int> random_order(n);
	std::iota(random_order.begin(), random_order.end(), 0);

	std::random_device rd;
	std::mt19937 gen(rd());

	do {
		changed_in_pass = false;

		std::shuffle(random_order.begin(), random_order.end(), gen);

		for (int i = 0; i < n; i++) {
			int v = random_order[i]; 

			bool moved = getBestDelta(g, v, partition);
			if (moved) {
				changed_in_pass = true;
				changed_global = true;
			}
		}
	} while (changed_in_pass);

	return changed_global;
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
	for (int i = 0; i < N; i++) {
		result[i] = i;
	}
	in_G = in;
	tot_G = tot;
}
template<class T>
louvain<T>::louvain(const graph<T>& G) {
	std::cout << std::setprecision(15);
	inizialization(G);
	bool flag = false;
	double old = getModularityOptimized(g, teck_partition, 0);
	first_modularity = old;
	std::cout << "Modularity: " << old << '\n';
	auto start = std::chrono::steady_clock::now();
	do {
		flag = moveNodes(g, teck_partition);
		if (flag) {
			aggregateGraph(g, teck_partition);
		}
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
	int n = g.getVertexCount();
	double m = g.getEdgeCount();
	double totalWeight = 1.0/(2.0 * m);
	community_count = teck_community_count;
	double result = 0;
	if (flag == 0) {
		for (int i = 0; i < community_count; i++) {
			result += in_G[i] * totalWeight - (tot_G[i] * totalWeight) * (tot_G[i] * totalWeight);
		}
	}
	else {
		for (int i = 0; i < community_count; i++) {
			result += in[i] * totalWeight - (tot[i] * totalWeight) * (tot[i] * totalWeight);
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
