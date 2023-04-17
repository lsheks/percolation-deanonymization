//
//  ethereum.h
//  ethereum
//
//  Created by Louis Shekhtman on 9/13/18.
//  Copyright © 2018 Louis Shekhtman. All rights reserved.
//

#include <stdio.h>
#include <map>
#include <unordered_set>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>

using std::vector;
using std::unordered_set;
using std::map;
using std::pair;
using std::string;
using std::unordered_map;

enum class BFS {Black, White, Gray};
class Ethereum {
  
private:
    uint N, num_links, S, crun, M, num_exchange_links;
    uint S_temp;
    unsigned long S_pairs;
    double k;
    vector< vector< uint > > adjacency_list;
    map< vector< uint >, bool > link_alive;
    vector< vector<uint> > all_links;
    vector< vector<uint> > all_links_share_noshare;
    vector< vector<uint> > all_links_share_share;
    vector< vector<uint> > all_links_noshare_share;

    vector< int > largest_component;
    vector < uint > component_sizes;
    vector < uint > nodes_order;
    vector < uint > distance_exchange;
    vector < uint > bfs_order_exchanges;
    vector<uint> exchange_nodes;
    vector<uint> link1;

    vector <bool> alive;
    vector <int> components;
    map <uint,uint> component_size;
    std::queue<uint> bfs_queue;
    vector <BFS> bfs_visited;
    
    std::mt19937 gen;
    
    std::uniform_int_distribution<uint> randint;
    std::uniform_int_distribution<uint> randint1;
    std::uniform_int_distribution<uint> randint2;
    vector< vector <double> > numlinks_S_history;
    vector< vector <double> > numlinks_S_history_modules;
   
public:
    int path_overall;
    Ethereum(uint r);
    void read_nodes();
    void read_nodes_dacunha();
    void read_nodes_dacunha_directed();
    void read_nodes_matjaz();
    void er_nodes(uint N, double k);
    void sf_nodes();
    void sf_nodes_clustering();
    void er_nodes_clustering();
    void read_nodes_sf_clustering();
    
    void order_by_deg();
    void link_exchange_last();
    void link_order_matjaz_dacunha(vector<uint> exchange_nodes);
    void link_exchange_er_last();
    void update_S( double param1);
    void update_S_directed( double param1);
    void find_S(double dominance);
    void BFS1();
    void BFS1_directed();

    void shuffle_links();
    void writeDist(uint runs, double perc_val, string text);
    void writeDist_directed(uint runs, double perc_val, double perc_val2, double perc_val3, string text);

    void perc_nodes(double perc_new);
    void perc_links(double perc_new);
    void perc_links_special(double perc_new1,double perc_new2,double perc_new3);

    void increase_p(double p_new);
    void writeS(uint runs, string text);
    void reach_node(uint node);
    void add_link(int sid, int tid);
    
    vector<vector< uint> > get_adjacency_list();
    vector<int> get_path_lengths(uint node);
    void read_er(vector <vector<uint> > adjacency_list1);
    
    vector<int> pathlength;
    
};
