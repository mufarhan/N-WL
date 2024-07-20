#ifndef NWL_WEISFEILERLEHMAN_H
#define NWL_WEISFEILERLEHMAN_H

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>

#include "Graph.h"

class NWL_WeisfeilerLehman {
  public:
    NWL_WeisfeilerLehman(string graph_database_name, int par_t, int par_d);
    NWL_WeisfeilerLehman();
    ~NWL_WeisfeilerLehman();

  private:
    // Compute lables for each graph in graph database.
    map<ulong, int> NeighborhoodWL(const Graph &g);
    void LoadIsomorphicTypes();
    vector<pair<int, int> > bfs_pair(int i, const Graph &g);

    // Reading a graph database from txt file.
    GraphDatabase read_graph_txt_file(string data_file);

    // Pairing function to map to a pair of Labels to a single label.
    ulong pairing(const ulong a, const ulong b);
    vector<int> g2, g3, g4, g5, g6, g7;
    int t, d;
};

#endif //KWL_WEISFEILERLEHMAN_H
