#include "louvain.h"
#include "dynamic_graph.h"
#include <string>

int main() {
	/*for (int i = 9; i < 28; i++) {
		std::string s = "../data/graph" + std::to_string(i) + ".txt";
		graph<std::vector<std::pair<int,double>>> G(s.data());

		std::cout << "Count of nodes: " << G.getVertexCount() << std::endl;
		std::cout << "Weight of edges = " << G.getEdgeCount() << std::endl;
		louvain<std::vector<std::pair<int, double>>> L(G);
		std::cout << std::endl;
	}*/

	//DynamicGraph G;
	//G.loadFromFile("../data/dynamic/graph2.txt");							
	//long long minTime = G.getMinTime(), maxTime = G.getMaxTime();
	//int k = 10;
	//long long period = (maxTime - minTime) / k;
	//for (long long i = minTime; i <= maxTime; i += period) {
	//	std::cout << "Period = " << i << ":" << i + period << std::endl;
	//	graph<std::vector<std::pair<int, double>>> snap = G.getSnapshot(i,i+period);
	//	//snap.printAdjList();
	//	std::cout << "Count of nodes: " << snap.getVertexCount() << std::endl;
	//	std::cout << "Weight of edges = " << snap.getEdgeCount() << std::endl;
	//	louvain<std::vector<std::pair<int, double>>> L(snap);
	//}

	std::string s = "../data/web-NotreDame.txt";
	graph<std::vector<std::pair<int, double>>> G(s.data());

	std::cout << "Count of nodes: " << G.getVertexCount() << std::endl;
	std::cout << "Weight of edges = " << G.getEdgeCount() << std::endl;
	louvain<std::vector<std::pair<int, double>>> L(G);
	std::cout << std::endl;
	return 0;
}