#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <utility>
#include<unordered_map>

template<class T>
class graph{
	int N;
	double M;
	std::vector<T> adj;
	std::vector<double> loops;
public:
	graph(std::vector <T> v = std::vector <T>(), std::vector<double> l = std::vector<double>(), int n = 0, int m = 0);
	graph(const char* filename, int flag = 0);
	graph(char* filename, char type);
	T & operator[](const int& i);
	const T& operator[](const int& i) const;
	double getWeightOfLoop(const int& v) const;
	int getVertexCount() const;
	int getEdgeCount() const;
	void printAdjList() const;
	void printNeighbours(const int& v) const;
};
template<class T>
graph<T>::graph(std::vector<T> v, std::vector <double> l, int n, int m) :N(n), M(m), adj(v), loops(l) {}
template<class T>
graph<T>::graph(const char* filename, int flag) { //если flag = 0, то в файле нет весов ребер, иначе  - есть
	std::ifstream file(filename, std::ios::in);
	file >> N >> M;
	int x, y;
	adj = std::vector<std::vector<std::pair<int, double>>>(N);
	loops = std::vector<double>(N, 0.0);
	std::vector<std::unordered_map<int, double>> adj_tmp = std::vector<std::unordered_map<int, double>>(N);
	if (flag == 0) {
		int count = 0;
		while (file >> x >> y) {
			if (x >= N || y >= N) {
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_int_distribution<> distrib(0, N - 1);

				x = distrib(gen);
				y = distrib(gen);
				while (x == y) { 
					y = distrib(gen);
				}
			}
			if (x != y) {
				auto result1 = adj_tmp[x].insert({ y,1.0 });
				auto result2 = adj_tmp[y].insert({ x,1.0 });
				if (result1.second) count++;
				if (result2.second) count++;
			}
		}
		//std::cout << count << std::endl;
		for (int i = 0; i < N; i++) { //добавляем петли сами с 0-ым весом
			adj_tmp[i].insert({ i,0.0 });
		}
		double m = 0;
		for (int v = 0; v < N; v++) {
			for (auto it = adj_tmp[v].begin(); it != adj_tmp[v].end(); it++) {
				adj[v].push_back({ (*it).first,(*it).second });
				m += (*it).second;
			}
		}
		//std::cout << m << std::endl;
		M = m / 2.0;
	}
}
template<class T>
graph<T>::graph(char* filename, char type) {
	std::ifstream finput;
	finput.open(filename, std::fstream::in);
	int x, y;
	finput >> x >> y;
	double nb_links = 0ULL;

	while (!finput.eof()) {
		unsigned int src, dest;
		double weight = 1.0L;

		finput >> src >> dest;
		
		if (finput) {
			if (adj.size() <= std::max(src, dest) + 1) {
				adj.resize(std::max(src, dest) + 1);
			}

			//adj[src].push_back(std::make_pair(dest, weight));
			if (src != dest) {
				adj[src].push_back(std::make_pair(dest, weight));
				adj[dest].push_back(std::make_pair(src, weight));
			}

			nb_links += 1ULL;
		}
	}

	finput.close();
	N = adj.size();
	M = nb_links;
	loops = std::vector<double>(N, 0.0);
}

template<class T>
T& graph<T>::operator[](const int& i) {
	return adj[i];
}
template<class T>
const T& graph<T>::operator[](const int& i) const {
	return adj[i];
}
template<class T>
int graph<T>::getVertexCount() const { return N; }
template<class T>
int graph<T>::getEdgeCount() const { return M; }
template<class T>
void graph<T>::printAdjList() const {
	for (int v = 0; v < adj.size(); v++) {
		std::cout << v << ": ";
		for (auto it = adj[v].begin(); it != adj[v].end(); it++) {
			std::cout << "{" << (*it).first << "," << (*it).second << "} ";
		}
		std::cout << std::endl;
	}
}
template<class T>
void graph<T>::printNeighbours(const int& v) const {
	std::cout << v << ": ";
	for (auto it = adj[v].begin(); it != adj[v].end(); it++) {
		std::cout << (*it).first << " ";
	}
	std::cout << std::endl;
}
template<class T>
double graph<T>::getWeightOfLoop(const int& v) const {
	return loops[v];
}
