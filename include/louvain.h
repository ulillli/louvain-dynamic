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
	std::vector<int> result;

	std::vector<double> neigh_weight; // neigh_weight[C] = dw(v, C) для текущей вершины v (-1 = C не среди соседей)
	std::vector<int> neigh_pos;       // сообщества, у вершин, являющихся соседями v
	int neigh_last = 0;               // кол-во соседей из разных сообществ

	graph<T> g;
	int teck_community_count = 0;
	std::vector<double> d;
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
	void computeNeighComm(const graph<T>& g, int v, const std::vector<int>& partition);
	void inizialization(const graph<T>& G);
	void reculculate(const graph<T>& g, const std::vector<int>& partition);
	void setSinglePartition(int n);
	double getModularity(const graph<T>& g, const std::vector<int>& partition);
	double getModularityOptimized(const graph<T>& g, const std::vector<int>& partition, bool flag);
	double getDw(int v, int C) const;
	double getGain(const graph<T>& g, const int& v, const std::vector<int>& partition, const int& C);
	void remove(int v, int C, const graph<T>& g, std::vector<int>& partition);
	void insert(int v, int C, const graph<T>& g, std::vector<int>& partition);
	bool getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition);
	bool aggregateGraph(graph<T>& g, std::vector<int>& partition);
	bool moveNodes(graph<T>& g, std::vector<int>& partition);
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
void louvain<T>::computeNeighComm(const graph<T>& g, int v, const std::vector<int>& partition) {
	for (int i = 0; i < neigh_last; i++) neigh_weight[neigh_pos[i]] = -1.0; // сбрасываем только то, что было измененно для проошлой вершины
	neigh_last = 0;
	for (const auto& e : g[v]) { 
		int u = e.first;
		if (u == v) continue;          
		int c = partition[u];
		if (neigh_weight[c] < 0) {     
			neigh_weight[c] = 0.0;
			neigh_pos[neigh_last++] = c;
		}
		neigh_weight[c] += e.second;
	}
}
template<class T>
void louvain<T>::reculculate(const graph<T>& g, const std::vector<int>& partition) {
	int n = g.getVertexCount();
	d = std::vector<double>(n, 0);
	neigh_weight.assign(n, -1.0);
	neigh_pos.assign(n, 0);
	neigh_last = 0;
	in = std::vector<double>(teck_community_count, 0);
	tot = std::vector<double>(teck_community_count, 0);
	for (int v = 0; v < n; v++) {
		int v_community = partition[v];
		in[v_community] += g.getWeightOfLoop(v);
		for (auto it = g[v].begin(); it != g[v].end(); it++) {
			int neighbour = (*it).first;
			int neighbour_community = partition[neighbour];
			tot[v_community] += (*it).second;
			d[v] += (*it).second;
			if (v != neighbour && neighbour_community == v_community) in[v_community] += (*it).second;
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
double louvain<T>::getDw(int v, int C) const {
	return neigh_weight[C] < 0 ? 0.0 : neigh_weight[C];
}
template<class T>
double louvain<T>::getGain(const graph<T>& g, const int& v, const std::vector<int>& partition, const int& C) {
	double totalWeight = 1.0 / (2.0 * g.getEdgeCount());
	double gain11 = (in[C] + 2.0 * getDw(v, C)) * totalWeight;
	double gain12 = ((tot[C] + d[v]) * totalWeight) * ((tot[C] + d[v]) * totalWeight);
	double gain21 = in[C] * totalWeight;
	double gain22 = (tot[C] * totalWeight) * (tot[C] * totalWeight);
	double gain23 = (d[v] * totalWeight) * (d[v] * totalWeight);
	return (gain11 - gain12) - (gain21 - gain22 - gain23);
}
template<class T>
void louvain<T>::remove(int v, int C, const graph<T>& g, std::vector<int>& partition) {
	in[C] = in[C] - 2.0 * getDw(v, C) - g.getWeightOfLoop(v);
	tot[C] = tot[C] - d[v];
	if (in[C] < 0) std::cout << "ERROR\n";
}
template<class T>
void louvain<T>::insert(int v, int C, const graph<T>& g, std::vector<int>& partition) {
	in[C] = in[C] + 2.0 * getDw(v, C) + g.getWeightOfLoop(v);
	tot[C] = tot[C] + d[v];
}
template<class T>
bool louvain<T>::getBestDelta(const graph<T>& g, const int& v, std::vector<int>& partition) {
	bool flag = false;
	int v_community = partition[v];
	computeNeighComm(g, v, partition);
	remove(v, v_community, g, partition);
	double best_gain = getGain(g, v, partition, v_community);
	int best_community = v_community;
	for (int i = 0; i < neigh_last; i++) {
		int community = neigh_pos[i];
		if (community != v_community) {
			double gain = getGain(g, v, partition, community);
			if (gain > best_gain) {
				best_gain = gain;
				best_community = community;
			}
		}
	}
	insert(v, best_community, g, partition);
	if (best_community != v_community) {
		flag = true;
		partition[v] = best_community;
	}
	return flag;
}
template<class T>
bool louvain<T>::aggregateGraph(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();

	std::vector<int> renumber(teck_community_count, -1);
	int new_count = 0;
	for (int v = 0; v < n; v++) {
		int c = partition[v];
		if (renumber[c] == -1) {
			renumber[c] = new_count++;
		}
	}
	if (new_count == n) return false; //������� ���: ���� �� ���������, partition/in/tot �������� ��������������
	for (int v = 0; v < n; v++) {
		partition[v] = renumber[partition[v]];
	}

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
	std::cout << "New graph with n = " << n << " and communities count = " << teck_community_count << std::endl;
	return true;
}
template<class T>
bool louvain<T>::moveNodes(graph<T>& g, std::vector<int>& partition) {
	int n = g.getVertexCount();
	std::cout << std::setprecision(15);
	bool flag = 1;
	bool changed = 0;
	double cur_modularity = getModularityOptimized(g, partition, 1);
	double new_modularity = cur_modularity;
	do {
		flag = 0;
		cur_modularity = new_modularity;
		for (int v = 0; v < n; v++) {
			if (getBestDelta(g, v, partition)) { //���� �� ���� ����������� �� ������ -> ����� ��� ������
				flag = true;
				changed = true;
			}
			/*if (best_delta.first > 0.0 && best_delta.second != v_community) {
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
			}*/
		}
		new_modularity = getModularityOptimized(g, partition, 1);
	} while (flag && new_modularity - cur_modularity > EPS); //����� ��������, ��� � ������ (�. 5.2)
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
	for (int i = 0; i < N; i++) {
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
		flag = moveNodes(g, teck_partition);
		if (flag) {
			flag = aggregateGraph(g, teck_partition);
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
	double m = g.getEdgeCount();
	community_count = teck_community_count;
	double result = 0;
	if (flag == 0) {
		for (int i = 0; i < (int)in_G.size(); i++) {
			result += in_G[i] / (2.0 * m) - (tot_G[i] / (2.0 * m)) * (tot_G[i] / (2.0 * m));
		}
	}
	else {
		for (int i = 0; i < community_count; i++) {
			result += in[i] / (2.0 * m) - (tot[i] / (2.0 * m)) * (tot[i] / (2.0 * m));
		}
	}
	return result;
}
template<class T>
double louvain<T>::getModularity(const graph<T>& g, const std::vector<int>& partition) {
	int community_count = 0;
	int n = g.getVertexCount();
	double m = g.getEdgeCount();
	for (int i = 0; i < n; i++) {
		if (partition[i] + 1 > community_count) community_count = partition[i] + 1;
	}
	std::vector<double> in = std::vector<double>(community_count, 0);
	std::vector<double> tot = std::vector<double>(community_count, 0);
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
		result += in[i] / (2.0 * m) - (tot[i] / (2.0 * m)) * (tot[i] / (2.0 * m));
	}
	return result;
}
