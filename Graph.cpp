#include "Graph.h"

  Graph::Graph() : m_num_nodes(0), m_num_edges(0) {}

  Graph::Graph(const uint num_nodes, const vector<pair<int, int> > &edgeList) {
    m_num_nodes = num_nodes;
    m_num_edges = edgeList.size();

    m_adjacency_lists.resize(num_nodes);

    for (auto const &e: edgeList) {
      add_edge(e.first, e.second);
    }
  }

  size_t Graph::add_node() {
    vector<int> new_node;
    m_adjacency_lists.push_back(move(new_node));
    m_num_nodes++;

    return m_num_nodes - 1;
  }

  void Graph::add_edge(const int v, const int w) {
    m_adjacency_lists[v].push_back(w);
    m_adjacency_lists[w].push_back(v);

    m_num_edges++;
  }

  void Graph::remove_edge(const int v, const int w) {
    m_adjacency_lists[v].erase(std::remove(m_adjacency_lists[v].begin(), m_adjacency_lists[v].end(), w), m_adjacency_lists[v].end());
    m_adjacency_lists[w].erase(std::remove(m_adjacency_lists[w].begin(), m_adjacency_lists[w].end(), v), m_adjacency_lists[w].end());

    m_num_edges--;
  }

  size_t Graph::get_degree(const int v) const {
    return m_adjacency_lists[v].size();
  }

  vector<int> Graph::get_neighbours(const int v) const {
    return m_adjacency_lists[v];
  }

  size_t Graph::get_num_nodes() const {
    return m_num_nodes;
  }

  size_t Graph::get_num_edges() const {
    return m_num_edges;
  }

  uint Graph::has_edge(const int v, const int w) const {
    if (find(m_adjacency_lists[v].begin(), m_adjacency_lists[v].end(), w) != m_adjacency_lists[v].end()) {
      return 1;
    } else {
      return 0;
    }
  }

  Graph::~Graph() {}
