#include "NWL_WeisfeilerLehman.h"

NWL_WeisfeilerLehman::NWL_WeisfeilerLehman() { }

NWL_WeisfeilerLehman::NWL_WeisfeilerLehman(string graph_database_name, int par_t, int par_d) {

      t = par_t;
      d = par_d;
      GraphDatabase gdb = read_graph_txt_file(graph_database_name);

      ulong num_graphs = gdb.size();
      vector<map<ulong, int> > color_counters;
      color_counters.reserve(num_graphs);
      //cout << num_graphs << endl;

      LoadIsomorphicTypes();
      for (auto &graph: gdb)
        color_counters.push_back(NeighborhoodWL(graph));

      int sum = 0;
      unordered_map<string, vector<ulong> > iso_test;
      for (ulong i = 0; i < color_counters.size(); i++) {
        map<ulong, int> c = color_counters[i];

        string hash_key = "";
        for (const auto &j: c) {
          ulong key = j.first;
          uint value = j.second;

          hash_key += to_string(key) + string("-") + to_string(value) + " ";
        }

        hash_key.pop_back();
        iso_test[hash_key].push_back(i);
      }

      std::unordered_map<string, vector<ulong> >::iterator it = iso_test.begin();
      while (it != iso_test.end()) {
        string key = it->first;
        vector<ulong> value = it->second;

        if(value.size() > 1)
          sum += (value.size() * (value.size() - 1)) / 2;
        it++;
      }
      cout << "Isomorphic Pairs = " << sum << endl;
}

vector<pair<int, int> > NWL_WeisfeilerLehman::bfs_pair(int i, const Graph &g) {

  size_t num_nodes = g.get_num_nodes();
  vector<pair<int, int> > neighborhood_subgraph;
  neighborhood_subgraph.push_back(make_pair(i, 0));

  queue<int> que;
  vector<int> P; P.resize(num_nodes, 99);

  que.push(i); P[i] = 0;
  while (!que.empty()) {
    int u = que.front();
    que.pop();

    vector<int> neighbors = g.get_neighbours(u);
    for (int w : neighbors) {
      if (P[w] == 99 && P[u] < d) {
        P[w] = P[u] + 1;
        neighborhood_subgraph.push_back(make_pair(w, P[w]));
        que.push(w);
      }
    }
  }

  return neighborhood_subgraph;
}

void NWL_WeisfeilerLehman::LoadIsomorphicTypes() {
  ifstream ifs("graphlists/graphlist2.txt"); int i;
  while(ifs >> i)
    g2.push_back(i);
  ifs.close();

  ifs.open("graphlists/graphlist3.txt");
  while(ifs >> i)
    g3.push_back(i);
  ifs.close();

  ifs.open("graphlists/graphlist4.txt");
  while(ifs >> i)
    g4.push_back(i);
  ifs.close();

  ifs.open("graphlists/graphlist5.txt");
  while(ifs >> i)
    g5.push_back(i);
  ifs.close();

  ifs.open("graphlists/graphlist6.txt");
  while(ifs >> i)
    g6.push_back(i);
  ifs.close();

  ifs.open("graphlists/graphlist7.txt");
  while(ifs >> i)
    g7.push_back(i);
  ifs.close();
}

map<ulong, int> NWL_WeisfeilerLehman::NeighborhoodWL(const Graph &g) {
  size_t num_nodes = g.get_num_nodes();

  vector<ulong> coloring(num_nodes);
  vector<ulong> coloring_temp(num_nodes);
  map<ulong, int> color_map, color_map_temp;

  // Assign isomorphism type to each 2-element set.
  for (int i = 0; i < num_nodes; ++i) {
    coloring[i] = 1;
    color_map[coloring[i]]++;
  }

  while (color_map.size() != color_map_temp.size()) {

    color_map_temp = color_map; color_map.clear();
    for (int i = 0; i < num_nodes; i++) {
      int temp_t = t;
      vector<pair<int, int> > neighborhood_subgraph = bfs_pair(i, g);

      int N = neighborhood_subgraph.size();
      if(t > N) temp_t = N;
      map<string, ulong> colors;
      if(temp_t == 1) {
	
	for (int v0 = 0; v0 < N; v0++) {
          vector<int> v; v.push_back(neighborhood_subgraph[v0].second); v.push_back(coloring[neighborhood_subgraph[v0].first]);
          colors["0"+to_string(neighborhood_subgraph[v0].second)] += pairing(v[0], v[1]);
	}

      } else if(temp_t == 2) {
	
        for (int v0 = 0; v0 < N-1; v0++) {
          for (int v1 = v0 + 1; v1 < N; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);

            vector<int> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(g2[code1]);
            sort(v.begin(), v.end());

	    ulong val = pairing(pairing(v[0], v[1]), v[2]);
            v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second);
            sort(v.begin(), v.end());
            string str = to_string(g2[code1]) + ", " + to_string(v[0]) + to_string(v[1]);
            colors[str] += val;
          }
        }
      } else if(temp_t == 3) {
        for (int v0 = 0; v0 < N-2; v0++) {
          for (int v1 = v0 + 1; v1 < N-1; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);
            for (int v2 = v1 + 1; v2 < N; v2++) {
	      int code2 = 2*code1 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v2].first);
	      code2 = 2*code2 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v2].first);

              vector<int> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(coloring[neighborhood_subgraph[v2].first]); v.push_back(g3[code2]);
              sort(v.begin(), v.end());

	      ulong val = pairing(pairing(pairing(v[0], v[1]), v[2]), v[3]);
	      v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second); v.push_back(neighborhood_subgraph[v2].second);
	      sort(v.begin(), v.end());
	      string str = to_string(g3[code2]) + ", " + to_string(v[0]) + to_string(v[1]) + to_string(v[2]);
              colors[str] += val;
	    }
          }
        }
      } else if (temp_t == 4) {

        for (int v0 = 0; v0 < N-3; v0++) {
          for (int v1 = v0 + 1; v1 < N-2; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);
            for (int v2 = v1 + 1; v2 < N-1; v2++) {
              int code2 = 2*code1 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v2].first);
              code2 = 2*code2 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v2].first);   // the order is important
              for (int v3 = v2 + 1; v3 < N; v3++) {
                int code3 = 2*code2 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v3].first);

                vector<ulong> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(coloring[neighborhood_subgraph[v2].first]); v.push_back(coloring[neighborhood_subgraph[v3].first]); v.push_back(g4[code3]);
		sort(v.begin(), v.end());

		ulong val = pairing(pairing(pairing(pairing(v[0], v[1]), v[2]), v[3]), v[4]);
		v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second); v.push_back(neighborhood_subgraph[v2].second); v.push_back(neighborhood_subgraph[v3].second);
		sort(v.begin(), v.end());
                string str = to_string(g4[code3]) + ", " + to_string(v[0]) + to_string(v[1]) + to_string(v[2]) + to_string(v[3]);
                colors[str] += val;
              }
            }
          }
        }
      } else if (temp_t == 5) {

        for (int v0 = 0; v0 < N-4; v0++) {
          for (int v1 = v0 + 1; v1 < N-3; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);
            for (int v2 = v1 + 1; v2 < N-2; v2++) {
              int code2 = 2*code1 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v2].first);
              code2 = 2*code2 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v2].first);   // the order is important
              for (int v3 = v2 + 1; v3 < N-1; v3++) {
                int code3 = 2*code2 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v3].first);
                for (int v4 = v3 + 1; v4 < N; v4++) {
                  int code4 = 2*code3 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v4].first);

                  vector<ulong> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(coloring[neighborhood_subgraph[v2].first]); v.push_back(coloring[neighborhood_subgraph[v3].first]); v.push_back(coloring[neighborhood_subgraph[v4].first]); v.push_back(g5[code4]);
                  sort(v.begin(), v.end());

		  ulong val = pairing(pairing(pairing(pairing(pairing(v[0], v[1]), v[2]), v[3]), v[4]), v[5]);
		  v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second); v.push_back(neighborhood_subgraph[v2].second); v.push_back(neighborhood_subgraph[v3].second); v.push_back(neighborhood_subgraph[v4].second);
		  sort(v.begin(), v.end());
                  string str = to_string(g5[code4]) + ", " + to_string(v[0]) + to_string(v[1]) + to_string(v[2]) + to_string(v[3]) + to_string(v[4]);
                  colors[str] += val;
                }
              }
            }
          }
        }
      } else if (temp_t == 6) {

        for (int v0 = 0; v0 < N-5; v0++) {
          for (int v1 = v0 + 1; v1 < N-4; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);
            for (int v2 = v1 + 1; v2 < N-3; v2++) {
              int code2 = 2*code1 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v2].first);
              code2 = 2*code2 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v2].first);   // the order is important
              for (int v3 = v2 + 1; v3 < N-2; v3++) {
                int code3 = 2*code2 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v3].first);
                for (int v4 = v3 + 1; v4 < N-1; v4++) {
                  int code4 = 2*code3 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v4].first);
                  for (int v5 = v4 + 1; v5 < N; v5++) {
                    int code5 = 2*code4 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v4].first, neighborhood_subgraph[v5].first);

                    vector<ulong> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(coloring[neighborhood_subgraph[v2].first]); v.push_back(coloring[neighborhood_subgraph[v3].first]); v.push_back(coloring[neighborhood_subgraph[v4].first]); v.push_back(coloring[neighborhood_subgraph[v5].first]); v.push_back(g6[code5]);
                    sort(v.begin(), v.end());

		    ulong val = pairing(pairing(pairing(pairing(pairing(pairing(v[0], v[1]), v[2]), v[3]), v[4]), v[5]), v[6]);
                    v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second); v.push_back(neighborhood_subgraph[v2].second); v.push_back(neighborhood_subgraph[v3].second); v.push_back(neighborhood_subgraph[v4].second); v.push_back(neighborhood_subgraph[v5].second);
		    sort(v.begin(), v.end());
                    string str = to_string(g6[code5]) + ", " + to_string(v[0]) + to_string(v[1]) + to_string(v[2]) + to_string(v[3]) + to_string(v[4]) + to_string(v[5]);
                    colors[str] += val;
                  }
                }
              }
	    }
          }
        }
      } else if (temp_t == 7) {

        for (int v0 = 0; v0 < N-6; v0++) {
          for (int v1 = v0 + 1; v1 < N-5; v1++) {
            int code1 = g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v1].first);
            for (int v2 = v1 + 1; v2 < N-4; v2++) {
              int code2 = 2*code1 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v2].first);
              code2 = 2*code2 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v2].first);   // the order is important
              for (int v3 = v2 + 1; v3 < N-3; v3++) {
                int code3 = 2*code2 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v3].first);
                code3 = 2*code3 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v3].first);
                for (int v4 = v3 + 1; v4 < N-2; v4++) {
                  int code4 = 2*code3 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v4].first);
                  code4 = 2*code4 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v4].first);
                  for (int v5 = v4 + 1; v5 < N-1; v5++) {
                    int code5 = 2*code4 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v5].first);
                    code5 = 2*code5 + g.has_edge(neighborhood_subgraph[v4].first, neighborhood_subgraph[v5].first);
                    for (int v6 = v5 + 1; v6 < N; v6++) {
                      int code6 = 2*code5 + g.has_edge(neighborhood_subgraph[v0].first, neighborhood_subgraph[v6].first);
                      code6 = 2*code6 + g.has_edge(neighborhood_subgraph[v1].first, neighborhood_subgraph[v6].first);
                      code6 = 2*code6 + g.has_edge(neighborhood_subgraph[v2].first, neighborhood_subgraph[v6].first);
                      code6 = 2*code6 + g.has_edge(neighborhood_subgraph[v3].first, neighborhood_subgraph[v6].first);
                      code6 = 2*code6 + g.has_edge(neighborhood_subgraph[v4].first, neighborhood_subgraph[v6].first);
                      code6 = 2*code6 + g.has_edge(neighborhood_subgraph[v5].first, neighborhood_subgraph[v6].first);

                      vector<ulong> v; v.push_back(coloring[neighborhood_subgraph[v0].first]); v.push_back(coloring[neighborhood_subgraph[v1].first]); v.push_back(coloring[neighborhood_subgraph[v2].first]); v.push_back(coloring[neighborhood_subgraph[v3].first]); v.push_back(coloring[neighborhood_subgraph[v4].first]); v.push_back(coloring[neighborhood_subgraph[v5].first]); v.push_back(coloring[neighborhood_subgraph[v6].first]); v.push_back(g7[code6]);
                      sort(v.begin(), v.end());

		      ulong val = pairing(pairing(pairing(pairing(pairing(pairing(pairing(v[0], v[1]), v[2]), v[3]), v[4]), v[5]), v[6]), v[7]);
                      v.clear(); v.push_back(neighborhood_subgraph[v0].second); v.push_back(neighborhood_subgraph[v1].second); v.push_back(neighborhood_subgraph[v2].second); v.push_back(neighborhood_subgraph[v3].second); v.push_back(neighborhood_subgraph[v4].second); v.push_back(neighborhood_subgraph[v5].second); v.push_back(neighborhood_subgraph[v6].second);
		      sort(v.begin(), v.end());
                      string str = to_string(g7[code6]) + ", " + to_string(v[0]) + to_string(v[1]) + to_string(v[2]) + to_string(v[3]) + to_string(v[4]) + to_string(v[5]) + to_string(v[6]);
                      colors[str] += val;
                    }
                  }
                }
              }
            }
          }
        }
      }

      vector<ulong> abc;
      for (auto m: colors) {
        abc.push_back(m.second);
	//cout << m.first << " " << m.second << endl;
      }
      sort(abc.begin(), abc.end());

      ulong new_color = coloring[i];
      for (long const &c: abc) {
          new_color = pairing(new_color, c);
      }

      coloring_temp[i] = new_color;
    }

    for (int i = 0 ; i < num_nodes ; i++) {
        coloring[i] = coloring_temp[i];

        color_map[coloring[i]]++;
    }
  }

  return color_map;
}

GraphDatabase NWL_WeisfeilerLehman::read_graph_txt_file(string data_file) {

  // Insert edges for each graph.
  ifstream edge_file(data_file); GraphDatabase graph_database;
  if (edge_file.is_open()) {

    int n, e, a, b, j = 0;
    while (edge_file >> n >> e) {
      vector<pair<int, int> > edge_list;

      Graph new_graph(n, edge_list);
      graph_database.push_back(new_graph);

      for(int i = 0; i < e; i++) {
        edge_file >> a >> b;

        graph_database[j].add_edge(a, b);
      } j++;
    }
    edge_file.close();
  } else {
    printf("%s", "!!! Unable to open file !!!\n");
    exit(EXIT_FAILURE);
  }

  return graph_database;
}

ulong NWL_WeisfeilerLehman::pairing(const ulong a, const ulong b) {
  return a >= b ? a * a + a + b : a + b * b;
}

NWL_WeisfeilerLehman::~NWL_WeisfeilerLehman() { }
