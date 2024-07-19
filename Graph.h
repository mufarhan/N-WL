#ifndef KWL_GRAPH_H
#define KWL_GRAPH_H

#include <map>
#include <vector>
#include <algorithm>

using namespace std;

namespace std {
    namespace {
        // Code from boost: Reciprocal of the golden ratio helps spread entropy and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine: http://stackoverflow.com/questions/4948780 .
        template<class T>
        inline void hash_combine(std::size_t &seed, T const &v) {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        // Recursive template code derived from Matthieu M.
        template<class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl {
            static void apply(size_t &seed, Tuple const &tuple) {
                HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
        };

        template<class Tuple>
        struct HashValueImpl<Tuple, 0> {
            static void apply(size_t &seed, Tuple const &tuple) {
                hash_combine(seed, get<0>(tuple));
            }
        };
    }

    template<typename ... TT>
    struct hash<std::tuple<TT...>> {
        size_t
        operator()(std::tuple<TT...> const &tt) const {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }
    };
}

class Graph {
  public:
    Graph();

    Graph(const uint num_nodes, const vector<pair<int, int> > &edgeList);

    // Add a single node to the graph.
    size_t add_node();

    // Add a single edge to the graph.
    void add_edge(const int v, const int w);

    // Remove a single edge to the graph.
    void remove_edge(const int v, const int w);

    // Get degree of node v.
    size_t get_degree(const int v) const;

    // Get neighbors of node v.
    vector<int> get_neighbours(const int v) const;

    // Get number of nodes in graph.
    size_t get_num_nodes() const;

    // Get number of edges in graph.
    size_t get_num_edges() const;

    // Returns 1 if edge {u,w} exists, otherwise 0.
    uint has_edge(const int v, const int w) const;

    ~Graph();

  private:
    vector<vector<int> > m_adjacency_lists;

    // Manage number of nodes in graph.
    size_t m_num_nodes;
    // Manage number of edges in graph.
    size_t m_num_edges;
};

typedef vector<Graph> GraphDatabase;

#endif //KWL_GRAPH_H
