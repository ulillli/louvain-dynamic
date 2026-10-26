#include "louvain.h"
//#include "dynamic_graph.h"
#include <string>

#ifndef DATA_DIR
#define DATA_DIR "./data"
#endif

int main() try {
	/*for (int i = 9; i < 28; i++) {
		std::string s = DATA_DIR "/graph" + std::to_string(i) + ".txt";
		graph<std::vector<std::pair<int,double>>> G(s.data());

		std::cout << "Count of nodes: " << G.getVertexCount() << std::endl;
		std::cout << "Weight of edges = " << G.getEdgeCount() << std::endl;
		louvain<std::vector<std::pair<int, double>>> L(G);
		std::cout << std::endl;
	}*/

	//DynamicGraph G;
	//G.loadFromFile(DATA_DIR "/dynamic/graph2.txt");							
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

	/*std::string s = DATA_DIR "/web-NotreDame.txt";
	graph<std::vector<std::pair<int, double>>> G(s.data());

	std::cout << "Count of nodes: " << G.getVertexCount() << std::endl;
	std::cout << "Weight of edges = " << G.getEdgeCount() << std::endl;
	louvain<std::vector<std::pair<int, double>>> L(G);
	std::cout << std::endl;*/

	for (int i = 1; i < 28; i++) {
		std::string s = DATA_DIR "/graph" + std::to_string(i) + ".txt";
		std::cout << s << std::endl;
		graph<std::vector<std::pair<int, double>>> G(s.data());

		std::cout << "Count of nodes: " << G.getVertexCount() << std::endl;
		std::cout << "Weight of edges = " << G.getEdgeCount() << std::endl;
		louvain<std::vector<std::pair<int, double>>> L(G);
		std::cout << std::endl;
	}
	return 0;
}
catch (const std::exception& e) {
	std::cerr << e.what() << std::endl;
	return 1;
}
