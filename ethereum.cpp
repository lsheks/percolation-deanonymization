//
//  ethereum.cpp
//  ethereum
//
//  Created by Louis Shekhtman on 9/13/18.
//  Copyright © 2018 Louis Shekhtman. All rights reserved.
//
#include "ethereum.h"
#include <algorithm>
#include<sstream>
#include <assert.h>
#include <cstdlib>
#include<string.h>
#include<iterator>
#include <unordered_map>


Ethereum::Ethereum(uint r){
        gen.seed((uint) time(0)+r);


}

void Ethereum::er_nodes(uint N, double k){
    S=N;
    M=k*N/2; //number of links
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);
    
    std::map<char,uint> node_ids;
    std::uniform_int_distribution<uint> randnode(0,N-1);
    uint tries=0,s,t;
    num_links=0;
    while(num_links < M){
        //std::cout<<num_links<<" out of "<< num_links_needed<<"  "<<mm<<"    "<<ll<<"    "<<k_vals[ll]<<"\n";
        do{
            tries+=1;
            s = randnode(gen);
            t = randnode(gen);
        }while(s == t);
        if(std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()){
            adjacency_list[s].push_back(t);
            adjacency_list[t].push_back(s);
            num_links++;
            if(s<t){
                    vector<uint> link={s, t};
                    link_alive[link]=1;
                    all_links.push_back(link);
                }
            if(s>t){
                    vector<uint> link={t, s};
                    link_alive[link]=1;
                    all_links.push_back(link);
                }
            }
        }

    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}
void Ethereum::sf_nodes(){
    N=10000;
    S=N;
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(2*N);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);
    
    std::map<char,uint> node_ids;
    std::uniform_int_distribution<uint> randnode(0,N-1);
    uint tries=0,s,t;
    num_links=0;
    
    vector<uint> stublist;
    uint num_times;
    double r;
    double kmin=3, kmax=N/3;
    double gamma=3;
    
    for(uint ii=0; ii<N; ++ii){
        r = ((double) rand() / (RAND_MAX));
        //std::cout<<"kmax/kmin="<<r*(pow(kmax/kmin,1-gamma)-1)+1<<"\n";
        num_times=kmin*pow(r*(pow(kmax/kmin,1-gamma)-1)+1,1/(1-gamma));
        for(uint jj=0; jj<num_times;++jj){
            stublist.push_back(ii);
        }
    }
    M=0;
    std::random_shuffle(stublist.begin(), stublist.end());
    std::cout<<"stublist size="<<stublist.size()<<"\n";
    for(uint uu=0; uu<(int) stublist.size()/2; ++uu){
            s = stublist[2*uu];
            t = stublist[2*uu+1];
            if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                adjacency_list[s].push_back(t);
                adjacency_list[t].push_back(s);
                M+=1;
               if(s<t){
                       vector<uint> link={s, t};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
               if(s>t){
                       vector<uint> link={t, s};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
        }
    }

    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}

void Ethereum::er_nodes_clustering(){
    N=5000;
    S=N;
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(2*N);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);
    
    std::map<char,uint> node_ids;
    std::uniform_int_distribution<uint> randnode(0,N-1);
    uint tries=0;
    num_links=0;

    
    uint num_times;
    vector<uint> node_deg(N,0);
    
    std::mt19937 generator(time(NULL));
    std::poisson_distribution<int> distribution(10);
    
    for(uint ii=0; ii<N; ++ii){
        //std::cout<<"kmax/kmin="<<r*(pow(kmax/kmin,1-gamma)-1)+1<<"\n";
        num_times=distribution(generator);
        node_deg[ii]=num_times;
    }
    std::cout<<"deg0="<<node_deg[0]<<"\n";
    
    vector<uint> remaining_node_edge_deg(node_deg);

    vector<uint> possible_triangles(N,0); //how many triangles this node could be involved in based on its degree
    vector<uint> triangles_stublist;
    for(uint ii=0; ii<N;++ii){
        possible_triangles[ii]=std::round(node_deg[ii]/2);
        for(uint jj=0; jj<possible_triangles[ii]; ++jj){
            triangles_stublist.push_back(ii);
        }
    }
    
    
    //Choose number of triangles node will be involved in
    uint N_triangles=10000;
    std::random_shuffle(triangles_stublist.begin(), triangles_stublist.end());
    uint n1, n2,n3,s ,t;
    vector<vector<uint> > node_pairs(3, vector<uint>(2,0));
    std::cout<<"triangles stublist size="<<triangles_stublist.size()<<"\n";
    uint n_t=0;
    uint ii=0;
    while(n_t<N_triangles and 3*ii+2<triangles_stublist.size()){
        n1=triangles_stublist[3*ii];
        n2=triangles_stublist[3*ii+1];
        n3=triangles_stublist[3*ii+2];
        ii+=1;
        if(n1!=n2 and n1!=n3 and n2!=n3){
            n_t+=1;
            node_pairs[0][0]=n1;
            node_pairs[0][1]=n2;
            node_pairs[1][0]=n1;
            node_pairs[1][1]=n3;
            node_pairs[2][0]=n2;
            node_pairs[2][1]=n3;
            for(auto pair :node_pairs){
                s=pair[0];
                t=pair[1];
                if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                    if (remaining_node_edge_deg[s]>0)
                        remaining_node_edge_deg[s]=remaining_node_edge_deg[s]-1;
                    if (remaining_node_edge_deg[s]==0)
                        std::cout<<"problem="<<s<<" "<<remaining_node_edge_deg[s]<<" "<<node_deg[s]<<"\n";
                    if (remaining_node_edge_deg[t]>0)
                        remaining_node_edge_deg[t]=remaining_node_edge_deg[t]-1;
                    
                    adjacency_list[s].push_back(t);
                    adjacency_list[t].push_back(s);
                    M+=1;
                   if(s<t){
                           vector<uint> link={s, t};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                   if(s>t){
                           vector<uint> link={t, s};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                    }
            }
        }
    }
    uint sum_of_degs = std::accumulate(node_deg.begin(), node_deg.end(), 0);
    uint sum_of_remaining = std::accumulate(remaining_node_edge_deg.begin(), remaining_node_edge_deg.end(), 0);
    std::cout<<"orig="<<sum_of_degs<<", remaining="<<sum_of_remaining<<"\n";
    vector<uint> stublist;
    for(uint ii=0; ii<N; ++ii){
        num_times=remaining_node_edge_deg[ii];
        for(uint jj=0; jj<num_times;++jj){
            stublist.push_back(ii);
        }
    }
    
    M=0;
    std::random_shuffle(stublist.begin(), stublist.end());
    std::cout<<"stublist size="<<stublist.size()<<"\n";
    for(uint uu=0; uu<(int) stublist.size()/2; ++uu){
            s = stublist[2*uu];
            t = stublist[2*uu+1];
            if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                adjacency_list[s].push_back(t);
                adjacency_list[t].push_back(s);
                M+=1;
               if(s<t){
                       vector<uint> link={s, t};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
               if(s>t){
                       vector<uint> link={t, s};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
        }
    }

    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    
    float clustering_coefficient, triplets=0,kval, triangles=0;

    vector<uint> neighbors;
    for(uint ii=0; ii<N;++ii){
        kval=adjacency_list[ii].size();
        triplets+=kval*(kval-1);
        std::sort(adjacency_list[ii].begin(), adjacency_list[ii].end());
        neighbors=adjacency_list[ii];
        for(uint jj=0; jj<neighbors.size(); ++jj){
            for (uint kk=0; kk<jj; ++kk){
                n1=neighbors[kk];
                n2=neighbors[jj];
                if(std::find(adjacency_list[n1].begin(), adjacency_list[n1].end(), n2) != adjacency_list[n1].end()){
                    triangles+=1;
                }
            }
        }
    }
    std::cout<<"triangles="<<triangles<<", clustering="<<triangles/triplets<<"\n";
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}
    
void Ethereum::sf_nodes_clustering(){
    N=10000;
    S=N;
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(2*N);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);
    
    std::map<char,uint> node_ids;
    std::uniform_int_distribution<uint> randnode(0,N-1);
    uint tries=0;
    num_links=0;

    
    uint num_times;
    double r;
    double kmin=3, kmax=N/4;
    double gamma=2.5;
    vector<uint> node_deg(N,0);
    for(uint ii=0; ii<N; ++ii){
        r = ((double) rand() / (RAND_MAX));
        //std::cout<<"kmax/kmin="<<r*(pow(kmax/kmin,1-gamma)-1)+1<<"\n";
        num_times=kmin*pow(r*(pow(kmax/kmin,1-gamma)-1)+1,1/(1-gamma));
        node_deg[ii]=num_times;
    }
    
    vector<uint> remaining_node_edge_deg(node_deg);

    vector<uint> possible_triangles(N,0); //how many triangles this node could be involved in based on its degree
    vector<uint> triangles_stublist;
    for(uint ii=0; ii<N;++ii){
        possible_triangles[ii]=std::round(node_deg[ii]/2);
        for(uint jj=0; jj<possible_triangles[ii]; ++jj){
            triangles_stublist.push_back(ii);
        }
    }
    
    
    //Choose number of triangles node will be involved in
    uint N_triangles=15000;
    std::random_shuffle(triangles_stublist.begin(), triangles_stublist.end());
    uint n1, n2,n3,s ,t;
    vector<vector<uint> > node_pairs(3, vector<uint>(2,0));
    std::cout<<"triangles stublist size="<<triangles_stublist.size()<<"\n";
    uint n_t=0;
    uint ii=0;
    while(n_t<N_triangles and 3*ii+2<triangles_stublist.size()){
        n1=triangles_stublist[3*ii];
        n2=triangles_stublist[3*ii+1];
        n3=triangles_stublist[3*ii+2];
        ii+=1;
        if(n1!=n2 and n1!=n3 and n2!=n3){
            n_t+=1;
            node_pairs[0][0]=n1;
            node_pairs[0][1]=n2;
            node_pairs[1][0]=n1;
            node_pairs[1][1]=n3;
            node_pairs[2][0]=n2;
            node_pairs[2][1]=n3;
            for(auto pair :node_pairs){
                s=pair[0];
                t=pair[1];
                if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                    if (remaining_node_edge_deg[s]==0)
                        std::cout<<"problem="<<s<<" "<<remaining_node_edge_deg[s]<<" "<<node_deg[s]<<" "<<std::count(triangles_stublist.begin(), triangles_stublist.end(), s)<<"\n";
                    if (remaining_node_edge_deg[s]>0)
                        remaining_node_edge_deg[s]=remaining_node_edge_deg[s]-1;
                    if (remaining_node_edge_deg[t]>0)
                        remaining_node_edge_deg[t]=remaining_node_edge_deg[t]-1;
                    
                    adjacency_list[s].push_back(t);
                    adjacency_list[t].push_back(s);
                    M+=1;
                   if(s<t){
                           vector<uint> link={s, t};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                   if(s>t){
                           vector<uint> link={t, s};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                    }
            }
        }
    }
    uint sum_of_degs = std::accumulate(node_deg.begin(), node_deg.end(), 0);
    uint sum_of_remaining = std::accumulate(remaining_node_edge_deg.begin(), remaining_node_edge_deg.end(), 0);
    std::cout<<"orig="<<sum_of_degs<<", remaining="<<sum_of_remaining<<"\n";
    vector<uint> stublist;
    for(uint ii=0; ii<N; ++ii){
        num_times=remaining_node_edge_deg[ii];
        for(uint jj=0; jj<num_times;++jj){
            stublist.push_back(ii);
        }
    }
    
    M=0;
    std::random_shuffle(stublist.begin(), stublist.end());
    std::cout<<"stublist size="<<stublist.size()<<"\n";
    for(uint uu=0; uu<(int) stublist.size()/2; ++uu){
            s = stublist[2*uu];
            t = stublist[2*uu+1];
            if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                adjacency_list[s].push_back(t);
                adjacency_list[t].push_back(s);
                M+=1;
               if(s<t){
                       vector<uint> link={s, t};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
               if(s>t){
                       vector<uint> link={t, s};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
        }
    }

    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    
    float clustering_coefficient, triplets=0,kval, triangles=0;

    vector<uint> neighbors;
    for(uint ii=0; ii<N;++ii){
        kval=adjacency_list[ii].size();
        triplets+=kval*(kval-1);
        std::sort(adjacency_list[ii].begin(), adjacency_list[ii].end());
        neighbors=adjacency_list[ii];
        for(uint jj=0; jj<neighbors.size(); ++jj){
            for (uint kk=0; kk<jj; ++kk){
                n1=neighbors[kk];
                n2=neighbors[jj];
                if(std::find(adjacency_list[n1].begin(), adjacency_list[n1].end(), n2) != adjacency_list[n1].end()){
                    triangles+=1;
                }
            }
        }
    }
    std::cout<<"triangles="<<triangles<<", clustering="<<triangles/triplets<<"\n";
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}
    


void Ethereum::read_nodes(){
    N=2721080;
    S=N;
    M=5262468; //number of links
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    for(uint i=0; i<2; ++i){ //change to i<10 for actually running
        std::ifstream inFile;
        std::ostringstream oss1;
        string line, node1, node2, delimiter=",";
        
        size_t pos;
        if(i==0){
            oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/mfg-trans.csv";
            //Change above path when running on cluster
            }
        if(i==1){
            oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/mfg-internal-trans.csv";
            }
        //std::cout<<(oss1.str()).c_str()<<"\n";
        inFile.open((oss1.str()).c_str());
        uint linenum=0;
        if (inFile.is_open()) {
            getline (inFile,line);
            //std::cout<<"line="<<line<<"\n";
            while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
                linenum+=1;
                //std::cout<<"line="<<line<<"\n";
                //std::cout<<linenum<<"   "<<i<<"\n";
                pos = line.find(delimiter);
                node1 = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                pos = line.find(delimiter);
                node2=line.substr(0,pos);
                
                uint converted1 = std::stoi(node1, nullptr, 10);
                uint converted2=std::stoi(node2, nullptr, 10);
                adjacency_list[converted1].push_back(converted2);
                adjacency_list[converted2].push_back(converted1);
                if(converted1<converted2){
                    vector<uint> link={converted1, converted2};
                    link_alive[link]=1;
                    all_links.push_back(link);
                }
                if(converted1>converted2){
                    vector<uint> link={converted2, converted1};
                    link_alive[link]=1;
                    all_links.push_back(link);

                }
            }
        }
        inFile.close();

    }
    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}


void Ethereum::read_er(vector <vector<uint> > adjacency_list1){
    N=adjacency_list1.size();
    S=N;
    M=0; //number of links
    for(uint aa=0; aa<adjacency_list1.size(); ++aa){
        M+=adjacency_list1[aa].size();
    }
    M=M/2;
    
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    uint node1, node2;
    for(uint i=0; i<adjacency_list1.size(); ++i){ //change to i<10 for actually running
        for(uint j=0; j<adjacency_list1[i].size(); ++j){
            node1=i;
            node2=adjacency_list1[i][j];
            if(node1<node2){
                adjacency_list[node1].push_back(node2);
                adjacency_list[node2].push_back(node1);
                vector<uint> link={node1, node2};
                link_alive[link]=1;
                all_links.push_back(link);
            }
        }
    }

    
    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}

void Ethereum::read_nodes_sf_clustering(){
    N=10000;
    S=N;
    M=820273; //number of links
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    
    vector<uint> triangles_stublist;
    std::ifstream inFile;
    std::ostringstream oss1;
    string line, node1, delimiter=",";
    size_t pos;
        oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/triangles_stublist081220.txt";
    inFile.open((oss1.str()).c_str());
    uint linenum=0;
    if (inFile.is_open()) {
        getline (inFile,line);
        while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
            linenum+=1;
            pos = line.find(delimiter);
            node1 = line.substr(0, pos);
            uint converted1 = std::stoi(node1, nullptr, 10);
            triangles_stublist.push_back(converted1);
        }
    }
    inFile.close();
    
    vector<uint> stublist;
    std::ifstream inFile1;
    std::ostringstream oss2;
        oss2<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/deg_stublist081220.txt";
    inFile1.open((oss2.str()).c_str());
    linenum=0;
    if (inFile1.is_open()) {
        getline (inFile1,line);
        while ( getline (inFile1,line)){ //remove and linnum<10000 when running for real
            linenum+=1;
            pos = line.find(delimiter);
            node1 = line.substr(0, pos);
            uint converted1 = std::stoi(node1, nullptr, 10);
            stublist.push_back(converted1);
        }
    }
    inFile1.close();

    

    
    //Choose number of triangles node will be involved in
    std::random_shuffle(triangles_stublist.begin(), triangles_stublist.end());
    uint n1, n2,n3,s ,t;
    vector<vector<uint> > node_pairs(3, vector<uint>(2,0));
    std::cout<<"triangles stublist size="<<triangles_stublist.size()<<"\n";
    uint n_t=0;
    uint ii=0;
    while(3*ii+2<triangles_stublist.size()){
        n1=triangles_stublist[3*ii];
        n2=triangles_stublist[3*ii+1];
        n3=triangles_stublist[3*ii+2];
        ii+=1;
        if(n1!=n2 and n1!=n3 and n2!=n3){
            n_t+=1;
            node_pairs[0][0]=n1;
            node_pairs[0][1]=n2;
            node_pairs[1][0]=n1;
            node_pairs[1][1]=n3;
            node_pairs[2][0]=n2;
            node_pairs[2][1]=n3;
            for(auto pair :node_pairs){
                s=pair[0];
                t=pair[1];
                if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                    adjacency_list[s].push_back(t);
                    adjacency_list[t].push_back(s);
                    M+=1;
                   if(s<t){
                           vector<uint> link={s, t};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                   if(s>t){
                           vector<uint> link={t, s};
                           link_alive[link]=1;
                           all_links.push_back(link);
                       }
                    }
            }
        }
    }
    std::random_shuffle(stublist.begin(), stublist.end());
    std::cout<<"stublist size="<<stublist.size()<<"\n";
    for(uint uu=0; uu<(int) stublist.size()/2; ++uu){
            s = stublist[2*uu];
            t = stublist[2*uu+1];
            if((std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) and s!=t){
                adjacency_list[s].push_back(t);
                adjacency_list[t].push_back(s);
                M+=1;
               if(s<t){
                       vector<uint> link={s, t};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
               if(s>t){
                       vector<uint> link={t, s};
                       link_alive[link]=1;
                       all_links.push_back(link);
                   }
        }
    }

    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    
    float clustering_coefficient, triplets=0,kval, triangles=0;

    vector<uint> neighbors;
    for(uint ii=0; ii<N;++ii){
        kval=adjacency_list[ii].size();
        triplets+=kval*(kval-1);
        std::sort(adjacency_list[ii].begin(), adjacency_list[ii].end());
        neighbors=adjacency_list[ii];
        for(uint jj=0; jj<neighbors.size(); ++jj){
            for (uint kk=0; kk<jj; ++kk){
                n1=neighbors[kk];
                n2=neighbors[jj];
                if(std::find(adjacency_list[n1].begin(), adjacency_list[n1].end(), n2) != adjacency_list[n1].end()){
                    triangles+=1;
                }
            }
        }
    }
    std::cout<<"triangles="<<triangles<<", clustering="<<triangles/triplets<<"\n";
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;
    
    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}

void Ethereum::read_nodes_dacunha_directed(){
    
    vector<uint> sharer_nodes{0, 1, 2, 3, 4, 4100, 10245, 7, 2049, 10249, 10, 2059, 4106, 13, 19, 21, 10266, 27, 4122, 29, 30, 8221, 2075, 2083, 4132, 37, 6181, 39, 4135, 10275, 6186, 43, 8236, 2094, 10287, 2097, 50, 51, 4145, 6197, 4153, 58, 59, 2106, 63, 4159, 8257, 64, 2115, 69, 6215, 2121, 10315, 76, 2123, 2126, 79, 10319, 4173, 8273, 10327, 88, 2136, 2138, 93, 2142, 95, 96, 8289, 6238, 10335, 2149, 6245, 8296, 2153, 2155, 108, 8300, 8303, 115, 2164, 6259, 4215, 4216, 2170, 6266, 2178, 2181, 138, 8330, 2187, 142, 8335, 8337, 10388, 8342, 4246, 156, 2205, 6302, 4257, 161, 2211, 4259, 4262, 6311, 2216, 8364, 179, 4278, 182, 2231, 4282, 4290, 4294, 200, 4296, 2250, 4302, 2260, 4309, 2267, 6364, 223, 2274, 6370, 8423, 6381, 6382, 2287, 2289, 8437, 4342, 248, 8443, 2302, 8457, 8458, 8459, 268, 2316, 270, 279, 2328, 8474, 6428, 4383, 8480, 291, 2339, 2343, 4391, 2344, 8490, 6447, 8497, 2354, 2360, 4410, 6460, 4422, 328, 8522, 4427, 2380, 2383, 2386, 8531, 342, 343, 8534, 2395, 2398, 6495, 357, 2405, 4457, 366, 8558, 2418, 2420, 8564, 2422, 4470, 2426, 8578, 2437, 4487, 395, 396, 4497, 2454, 2455, 411, 8607, 8617, 4523, 429, 430, 2479, 8624, 436, 4533, 438, 4535, 442, 445, 8637, 6591, 448, 8638, 2498, 8643, 451, 453, 454, 2503, 2504, 2501, 2510, 465, 2516, 469, 4568, 2522, 6620, 477, 8672, 8677, 489, 8682, 6633, 6634, 6641, 8690, 4595, 4597, 502, 8698, 508, 4609, 514, 4617, 532, 2586, 6684, 4639, 545, 4645, 8741, 8742, 8744, 6700, 6701, 8750, 2611, 6712, 2620, 6727, 8777, 590, 2640, 592, 6751, 6752, 4709, 8808, 8815, 4722, 8819, 4728, 6778, 637, 4735, 649, 657, 4761, 668, 8862, 4767, 8866, 4776, 6828, 8877, 2737, 4785, 2740, 698, 8891, 6845, 6847, 6849, 4802, 6855, 8904, 8907, 6861, 2768, 6864, 4818, 6870, 727, 4825, 4831, 737, 741, 2791, 748, 751, 755, 757, 6905, 763, 4871, 8973, 2832, 8978, 6936, 6939, 4897, 802, 803, 804, 4911, 9009, 2867, 4916, 6965, 822, 4919, 4920, 2874, 9025, 6977, 4929, 4932, 9028, 840, 4937, 9035, 4940, 2897, 2900, 864, 7008, 865, 7011, 4964, 7013, 9062, 2921, 2932, 2938, 904, 5001, 2956, 5005, 5004, 7057, 916, 921, 9113, 2975, 928, 7072, 2979, 9124, 9126, 5030, 2985, 7083, 7084, 9133, 2994, 7101, 7104, 5061, 7110, 967, 968, 7112, 970, 5068, 973, 974, 7123, 9172, 7129, 3033, 992, 7142, 9191, 7152, 3062, 9216, 5123, 1032, 1037, 7182, 9230, 7184, 5136, 1038, 7191, 1048, 7198, 7202, 7205, 5159, 1066, 3117, 3119, 7216, 9263, 9271, 9277, 7232, 7234, 1091, 9284, 9287, 9288, 1101, 1104, 3156, 7256, 9310, 1120, 5222, 9321, 3179, 7276, 1140, 3188, 3190, 9332, 5240, 1144, 9337, 1147, 9340, 1149, 1150, 7286, 1152, 5249, 5254, 3214, 3216, 5265, 1174, 1177, 3226, 7321, 5280, 3234, 5287, 5298, 5300, 7350, 1210, 1211, 1223, 7367, 1227, 1236, 9434, 1246, 1247, 7391, 5346, 3304, 5365, 7414, 5366, 1273, 1274, 9474, 7426, 7428, 7429, 5392, 7445, 7450, 1311, 9520, 1329, 1330, 5426, 3379, 9526, 9527, 5430, 3384, 3387, 3395, 5443, 3398, 7495, 1351, 5452, 1358, 3412, 3423, 5472, 9575, 9579, 5484, 9583, 5488, 5490, 5492, 1399, 3450, 5498, 1405, 3453, 1419, 1421, 3470, 1425, 7571, 5525, 1430, 7576, 7578, 3483, 1438, 9630, 1440, 1443, 1445, 5548, 1456, 3506, 1460, 3508, 1467, 3522, 9667, 1477, 3528, 1481, 9674, 1483, 7626, 1485, 3534, 3538, 3542, 5594, 5596, 7645, 9693, 9696, 1505, 3552, 3562, 5617, 3569, 1526, 1527, 5625, 9722, 1536, 3589, 9735, 9739, 5646, 1553, 9746, 5650, 3608, 7704, 3617, 7721, 7726, 3636, 5688, 7737, 9784, 5694, 5703, 9801, 5706, 3659, 7756, 3662, 5711, 3664, 9808, 9807, 5715, 9812, 3669, 7767, 7770, 5723, 9820, 5725, 7789, 1646, 1647, 1648, 1645, 1654, 1658, 1659, 1662, 1664, 3713, 1667, 9860, 1669, 1670, 1671, 1672, 3719, 1674, 3723, 3729, 1681, 1682, 1686, 9882, 3741, 1694, 1695, 1698, 3747, 5797, 1704, 1708, 1709, 9902, 7855, 3759, 5805, 1713, 7859, 1720, 3769, 3768, 1724, 1725, 1728, 7872, 1730, 1729, 1732, 1733, 1741, 3795, 3801, 1758, 7905, 5858, 1767, 3816, 7916, 5870, 1776, 7922, 3826, 5878, 7928, 1786, 9979, 1790, 9984, 1793, 3844, 5896, 1802, 7947, 5900, 5905, 5906, 3862, 3863, 7958, 10009, 1816, 3864, 10011, 1823, 3872, 5922, 5923, 3876, 3880, 5929, 10026, 1836, 5933, 3887, 1843, 5940, 5941, 1846, 7996, 10047, 1856, 3907, 1862, 1867, 1872, 3928, 1881, 1885, 8030, 1889, 3939, 1892, 3941, 1894, 1895, 10088, 10094, 3953, 8054, 6009, 3963, 8060, 8062, 3977, 8082, 3986, 1942, 8090, 1947, 3994, 4001, 8104, 4010, 10155, 4011, 10154, 1967, 4019, 4020, 1978, 1980, 8126, 1983, 10175, 10178, 1990, 6088, 1998, 4047, 1999, 6098, 4051, 8146, 10197, 10196, 8155, 2013, 4061, 6111, 10207, 2023, 10216, 2030, 2038, 2039, 2042, 6140};
    N=10407;
    S=N;
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links_share_noshare.reserve(M);
    all_links_noshare_share.reserve(M);
    all_links_share_share.reserve(M);

    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    for(uint i=0; i<1; ++i){ //change to i<10 for actually running
        std::ifstream inFile;
        std::ostringstream oss1;
        string line, node1, node2, delimiter=",";
        
        size_t pos;
        if(i==0){
            oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/data-da-cunha-full.csv";
            //Change above path when running on cluster
            }
      
        //std::cout<<(oss1.str()).c_str()<<"\n";
        inFile.open((oss1.str()).c_str());
        uint linenum=0;
        uint added=0;
        if (inFile.is_open()) {
            getline (inFile,line);
            //std::cout<<"line="<<line<<"\n";
            while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
                linenum+=1;
                //std::cout<<"line="<<line<<"\n";
                //std::cout<<linenum<<"   "<<i<<"\n";
                pos = line.find(delimiter);
                node1 = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                pos = line.find(delimiter);
                node2=line.substr(0,pos);
                
                uint converted1 = std::stoi(node1, nullptr, 10);
                uint converted2=std::stoi(node2, nullptr, 10);
                if((std::find(sharer_nodes.begin(), sharer_nodes.end(), converted1) == sharer_nodes.end())
                   && (std::find(adjacency_list[converted2].begin(),adjacency_list[converted2].end(), converted1)==adjacency_list[converted2].end())){
                    adjacency_list[converted2].push_back(converted1);
                    added+=1;
                    vector<uint> link={converted1, converted2};
                    link_alive[link]=1;
                    all_links_share_noshare.push_back(link);

                    adjacency_list[converted1].push_back(converted2);
                    added+=1;
                    vector<uint> link2={converted2, converted1};
                    link_alive[link2]=1;
                    all_links_noshare_share.push_back(link2);
                }
                if((std::find(sharer_nodes.begin(), sharer_nodes.end(), converted1) != sharer_nodes.end())
                   && (std::find(adjacency_list[converted1].begin(),adjacency_list[converted1].end(), converted2)==adjacency_list[converted1].end())){
                    adjacency_list[converted2].push_back(converted1);
                    added+=1;
                    vector<uint> link={converted1, converted2};
                    link_alive[link]=1;
                    all_links_share_share.push_back(link);

                    adjacency_list[converted1].push_back(converted2);
                    added+=1;
                    vector<uint> link2={converted2, converted1};
                    link_alive[link2]=1;
                    all_links_share_share.push_back(link2);
                }
            }
        }
        inFile.close();
        std::cout<<"added="<<added<<"\n";

    }
    M=all_links_share_share.size()+all_links_share_noshare.size()+all_links_noshare_share.size(); //number of links

    std::shuffle(all_links_share_share.begin(),all_links_share_share.end(), gen);
    std::shuffle(all_links_share_noshare.begin(),all_links_share_noshare.end(), gen);
    std::shuffle(all_links_noshare_share.begin(),all_links_noshare_share.end(), gen);

    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}


void Ethereum::read_nodes_dacunha(){
    N=10407;
    S=N;
    M=820273; //number of links
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    for(uint i=0; i<1; ++i){ //change to i<10 for actually running
        std::ifstream inFile;
        std::ostringstream oss1;
        string line, node1, node2, delimiter=",";
        
        size_t pos;
        if(i==0){
            oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/data-da-cunha-full.csv";
            //Change above path when running on cluster
            }
      
        //std::cout<<(oss1.str()).c_str()<<"\n";
        inFile.open((oss1.str()).c_str());
        uint linenum=0;
        uint added=0;
        if (inFile.is_open()) {
            getline (inFile,line);
            //std::cout<<"line="<<line<<"\n";
            while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
                linenum+=1;
                //std::cout<<"line="<<line<<"\n";
                //std::cout<<linenum<<"   "<<i<<"\n";
                pos = line.find(delimiter);
                node1 = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                pos = line.find(delimiter);
                node2=line.substr(0,pos);
                
                uint converted1 = std::stoi(node1, nullptr, 10);
                uint converted2=std::stoi(node2, nullptr, 10);
                if(std::find(adjacency_list[converted1].begin(), adjacency_list[converted1].end(), converted2) == adjacency_list[converted1].end()){
                    adjacency_list[converted1].push_back(converted2);
                    adjacency_list[converted2].push_back(converted1);
                    added+=1;
                    if(converted1<converted2){
                        vector<uint> link={converted1, converted2};
                        link_alive[link]=1;
                        all_links.push_back(link);
                    }
                    if(converted1>converted2){
                        vector<uint> link={converted2, converted1};
                        link_alive[link]=1;
                        all_links.push_back(link);
                    }
                }
            }
        }
        inFile.close();
        std::cout<<"added="<<added<<"\n";

    }
    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}

vector<vector< uint> > Ethereum::get_adjacency_list(){
    return (adjacency_list);
}

void Ethereum::read_nodes_matjaz(){
    N=309;
    S=N;
    M=3282; //number of links
    adjacency_list.resize(N); //based on number of nodes in dataset for ethereum
    alive.resize(N);
    all_links.reserve(M);
    components.resize(N);
    largest_component.resize(N);
    bfs_visited.resize(N);
    nodes_order.resize(N);
    distance_exchange=vector<uint>(N, 5000000);
    std::uniform_int_distribution<uint> randint(0,N-1);
    std::fill(alive.begin(),alive.end(),true);

    std::map<char,uint> node_ids;
    for(uint i=0; i<1; ++i){
        std::ifstream inFile;
        std::ostringstream oss1;
        string line, node1, node2, delimiter=",";
        
        size_t pos;
        if(i==0){
            oss1<<"/Users/louisshekhtman/Documents/ethereum/networkEdge/data-matjaz-corruption-070420.csv";
            //Change above path when running on cluster
            }
        //std::cout<<(oss1.str()).c_str()<<"\n";
        inFile.open((oss1.str()).c_str());
        uint linenum=0;
        if (inFile.is_open()) {
            getline (inFile,line);
            //std::cout<<"line="<<line<<"\n";
            while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
                linenum+=1;
                //std::cout<<"line="<<line<<"\n";
                //std::cout<<linenum<<"   "<<i<<"\n";
                pos = line.find(delimiter);
                node1 = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                pos = line.find(delimiter);
                node2=line.substr(0,pos);
                
                uint converted1 = std::stoi(node1, nullptr, 10);
                uint converted2=std::stoi(node2, nullptr, 10);
                adjacency_list[converted1].push_back(converted2);
                adjacency_list[converted2].push_back(converted1);
                if(converted1<converted2){
                    vector<uint> link={converted1, converted2};
                    link_alive[link]=1;
                    all_links.push_back(link);
                }
                if(converted1>converted2){
                    vector<uint> link={converted2, converted1};
                    link_alive[link]=1;
                    all_links.push_back(link);

                }
            }
        }
        inFile.close();

    }
    std::shuffle(all_links.begin(),all_links.end(), gen);
    //std::cout<<"all_links size="<<all_links.size()<<"\n";

    for(uint ii=0; ii<N; ++ii){
        nodes_order[ii]=ii;
    }
    std::shuffle(nodes_order.begin(),nodes_order.begin()+N, gen);
    //std::cout<<"nodes_order[99934]="<<nodes_order[99934]<<"\n";

    //if we aren't using the exchanges first then this doesn't matter,
    //but for when we want the distance this will be useful to have in order.
    bfs_order_exchanges=nodes_order;

}


void Ethereum::order_by_deg(){
        std::ifstream inFile;
        std::ostringstream oss1;
        string line, node1, node2, delimiter=",";
    
        size_t pos;
        oss1<<"highest_degree.csv";
        std::cout<<(oss1.str()).c_str()<<"\n";
        inFile.open((oss1.str()).c_str());
        uint linenum=0;
        nodes_order.clear();
        if (inFile.is_open()) {
            getline (inFile,line);
            //std::cout<<"line="<<line<<"\n";
            while ( getline (inFile,line)){ //remove and linnum<10000 when running for real
                linenum+=1;
                pos = line.find(delimiter);
                node1 = line.substr(0, pos);
                uint converted1 = std::stoi(node1, nullptr, 10);
                nodes_order.push_back(converted1);
            }
        }

        inFile.close();
    
        vector<bool> added(N, false);
        for(auto i : nodes_order){
            added[i]=true;
        }
        uint count=0;
        for(uint i=0; i<N; ++i){
            if(added[i]==false){
                    nodes_order.push_back(i);
                    count+=1;
                    std::cout<<count<<"\n";
            }
        }
}

void Ethereum::link_exchange_er_last(){
       //Here just choose random node to be exchange
    exchange_nodes={0}; //all exchanges
    std::uniform_int_distribution<uint> randint(0,N-1);
    uint enode=randint(gen);
    exchange_nodes[0]=enode;
    vector< vector<uint> > exchange_links;
    vector< vector<uint> > non_exchange_links;
    exchange_links.reserve(M);
    non_exchange_links.reserve(M);
    bool is_exchange_link;
    for ( auto ll : all_links){
        is_exchange_link=false;
        /*for ( auto nn : exchange_nodes){
            if (ll[0]==nn || ll[1]==nn){
                exchange_links.push_back(ll);
                is_exchange_link=true;
                break;
            }
        }*/
        if(is_exchange_link==false){
            non_exchange_links.push_back(ll);
        }
        
    }
    num_exchange_links=exchange_links.size();
    vector< vector<uint> > all_links2;
    all_links2.reserve(M);
    all_links2.insert(all_links2.end(), non_exchange_links.begin(), non_exchange_links.end());
    all_links2.insert(all_links2.end(), exchange_links.begin(), exchange_links.end());

    all_links=all_links2;
    
    //Carry out the BFS based using the exchanges first
    for(uint nn=0; nn<exchange_nodes.size(); ++nn){
        auto it_nn= find (bfs_order_exchanges.begin(), bfs_order_exchanges.end(), exchange_nodes[nn]);
        std::iter_swap(bfs_order_exchanges.begin()+nn , it_nn);
        std::cout<<"position "<<nn<<"="<<bfs_order_exchanges[nn]<<"\n";
        }

}

void Ethereum::shuffle_links(){
    //Build stublist
    uint check_node=304;
    std::cout<<check_node<<": ";
    for (auto nn : adjacency_list[check_node]){
        std::cout<<nn<<" ";
    }
    std::cout<<"\n";
    vector<uint> stub_list(2*M,0);
    std::cout<<"num_links="<<M<<"\n";
    uint count=0;
    for (auto link : all_links){
        stub_list[count]=link[0];
        stub_list[count+1]=link[1];
        count+=2;
    }
    std::cout<<"num_links="<<count<<" (should be twice above value, due to double count of stubs)\n";
    
    //Shuffle links
    std::shuffle(stub_list.begin(),stub_list.end(), gen);
    vector< vector< uint > > adjacency_list2;
    adjacency_list2.resize(N); //based on number of nodes in dataset for ethereum
    map< vector< uint >, bool > link_alive2;
    vector< vector<uint> > all_links2;
    all_links2.reserve(M);
    uint node1,node2;
    for(uint ii=0; ii<stub_list.size(); ii+=2){
        node1=stub_list[ii];
        node2=stub_list[ii+1];
        if(node1!=node2 and std::find(adjacency_list2[node1].begin(), adjacency_list2[node1].end(), node2) == adjacency_list2[node1].end()){
            adjacency_list2[node1].push_back(node2);
            adjacency_list2[node2].push_back(node1);
            if(node1<node2){
                vector<uint> link={node1, node2};
                link_alive2[link]=1;
                all_links2.push_back(link);
            }
            if(node1>node2){
                vector<uint> link={node2, node1};
                link_alive2[link]=1;
                all_links2.push_back(link);
            }
        }
    }
    adjacency_list=adjacency_list2;
    all_links=all_links2;
    link_alive=link_alive2;
    std::cout<<check_node<<": ";
    for (auto nn : adjacency_list[check_node]){
        std::cout<<nn<<" ";
    }
    std::cout<<"\n";
    
}


void Ethereum::link_exchange_last(){
       //Node ids for different exchanges:
       //{'shapeshift': [2319667, 1025299, 2472762, 1744523],
       //'poloniex': [1142452],
       //'kraken': [1156123], 'bittrex': [2098102], 'changelly': [2248226]}

        //exchange_nodes={2319667, 1025299, 2472762, 1744523}; //shapeshift
	//exchange_nodes={1142452}; //poloniex
	//exchange_nodes={1156123}; //kraken
        //exchange_nodes={2098102}; //bittrex
	//exchange_nodes={2248226}; //changelly
	exchange_nodes={2319667,1025299, 2472762,1744523,1142452,1156123,2098102,2248226}; //all exchanges	
	vector< vector<uint> > exchange_links;
        vector< vector<uint> > non_exchange_links;
        exchange_links.reserve(M);
        non_exchange_links.reserve(M);
        bool is_exchange_link;
        for ( auto ll : all_links){
            is_exchange_link=false;
            for ( auto nn : exchange_nodes){
                if (ll[0]==nn || ll[1]==nn){
                    exchange_links.push_back(ll);
                    is_exchange_link=true;
                    break;
                }
            }
            if(is_exchange_link==false){
                non_exchange_links.push_back(ll);
            }
            
        }
        num_exchange_links=exchange_links.size();
	vector< vector<uint> > all_links2;
        all_links2.reserve(M);
        all_links2.insert(all_links2.end(), non_exchange_links.begin(), non_exchange_links.end());
        all_links2.insert(all_links2.end(), exchange_links.begin(), exchange_links.end());

        all_links=all_links2;
    
        //Carry out the BFS based using the exchanges first
        for(uint nn=0; nn<exchange_nodes.size(); ++nn){
            auto it_nn= find (bfs_order_exchanges.begin(), bfs_order_exchanges.end(), exchange_nodes[nn]);
            std::iter_swap(bfs_order_exchanges.begin()+nn , it_nn);
            std::cout<<"position "<<nn<<"="<<bfs_order_exchanges[nn]<<"\n";
            }

}

void Ethereum::link_order_matjaz_dacunha(vector<uint> exchange_nodes1){

    exchange_nodes=exchange_nodes1; //all exchanges
    //std::cout<<"exchnage="<<exchange_nodes[0]<<"\n";
    vector< vector<uint> > exchange_links;
        vector< vector<uint> > non_exchange_links;
        exchange_links.reserve(M);
        non_exchange_links.reserve(M);
        bool is_exchange_link;
        for ( auto ll : all_links){
            is_exchange_link=false;
            if(is_exchange_link==false){
                non_exchange_links.push_back(ll);
            }
            
        }
    num_exchange_links=0;// for these networks there are no exchange links
    std::shuffle(non_exchange_links.begin(),non_exchange_links.end(), gen);
    //std::cout<<"len="<<non_exchange_links.size()<<"\n";

    vector< vector<uint> > all_links2;
        all_links2.reserve(M);
        all_links2.insert(all_links2.end(), non_exchange_links.begin(), non_exchange_links.end());

        all_links=all_links2;
    
        //Carry out the BFS based using the exchanges first
        for(uint nn=0; nn<exchange_nodes.size(); ++nn){
            auto it_nn= find (bfs_order_exchanges.begin(), bfs_order_exchanges.end(), exchange_nodes[nn]);
            std::iter_swap(bfs_order_exchanges.begin()+nn , it_nn);
            std::cout<<"position "<<nn<<"="<<bfs_order_exchanges[nn]<<"\n";
            }
    vector<uint> nodes_order(N,0);
    nodes_order[0]=exchange_nodes1[0];
    uint place_node=0;
    for(uint ii=1; ii<nodes_order.size();++ii){
        if(place_node==exchange_nodes1[0]){
            place_node+=1;
        }
        nodes_order[ii]=place_node;
        place_node+=1;
    }
    bfs_order_exchanges=nodes_order;

}


void Ethereum::add_link(int sid, int tid){
 // std::cout << "Adding link " <<sid << " -- " << tid << std::endl;
if(std::find(adjacency_list[sid].begin(), adjacency_list[sid].end(), tid) == adjacency_list[sid].end()){
	    adjacency_list[sid].push_back(tid);
	    adjacency_list[tid].push_back(sid);
	}
}

void Ethereum::update_S_directed(double param1)
{
    std::fill(distance_exchange.begin(), distance_exchange.end(),5000000);
    BFS1_directed();
    std::pair<uint, uint> LcbarIndexPair;
    auto LcbarIndexPair1=std::max_element(component_size.begin(), component_size.end(), [](std::pair<const uint,uint>& a, std::pair<const uint,uint>&b){return a.second < b.second;});
    LcbarIndexPair=std::pair<uint, uint>(LcbarIndexPair1->first, LcbarIndexPair1->second);

    uint giant_index = LcbarIndexPair.first;
    component_sizes.clear();
    for(map<uint,uint>::iterator it = component_size.begin(); it != component_size.end(); ++it) {
      component_sizes.push_back(it->second);
      //std::cout << it->second << "\n";
    }
    std::sort(component_sizes.rbegin(), component_sizes.rend());
    
    if(std::count(components.begin(), components.end(), giant_index) != LcbarIndexPair.second){
        std::cerr << "Miscounted components!! : counted " <<  LcbarIndexPair.second << " found " <<
        std::count(components.begin(), components.end(), giant_index) << "!\n";
    }
    S=LcbarIndexPair.second;
    
    std::vector<double> to_add;
    to_add.push_back(param1);
    for(uint ii=0; ii<100; ++ii){ //count 100 largest components
       if(ii<component_sizes.size()){
            to_add.push_back(component_sizes[ii]);
        }
    }
    numlinks_S_history.emplace_back(to_add);
}



void Ethereum::update_S(double param1)
{
    std::fill(distance_exchange.begin(), distance_exchange.end(),5000000);
    BFS1();
    std::pair<uint, uint> LcbarIndexPair;
    auto LcbarIndexPair1=std::max_element(component_size.begin(), component_size.end(), [](std::pair<const uint,uint>& a, std::pair<const uint,uint>&b){return a.second < b.second;});
    LcbarIndexPair=std::pair<uint, uint>(LcbarIndexPair1->first, LcbarIndexPair1->second);

    uint giant_index = LcbarIndexPair.first;
    component_sizes.clear();
    for(map<uint,uint>::iterator it = component_size.begin(); it != component_size.end(); ++it) {
      component_sizes.push_back(it->second);
      //std::cout << it->second << "\n";
    }
    std::sort(component_sizes.rbegin(), component_sizes.rend());
    
    if(std::count(components.begin(), components.end(), giant_index) != LcbarIndexPair.second){
        std::cerr << "Miscounted components!! : counted " <<  LcbarIndexPair.second << " found " <<
        std::count(components.begin(), components.end(), giant_index) << "!\n";
    }
    S=LcbarIndexPair.second;
    
    std::vector<double> to_add;
    to_add.push_back(param1);
    for(uint ii=0; ii<100; ++ii){ //count 100 largest components
       if(ii<component_sizes.size()){
            to_add.push_back(component_sizes[ii]);
        }
    }
    numlinks_S_history.emplace_back(to_add);
}


void Ethereum::perc_nodes(double perc_new){
    uint num_to_remove= perc_new*N;
    for (uint uu=0; uu<num_to_remove; ++uu){
        uint node_to_remove=nodes_order[uu];
        alive[node_to_remove]=false;

    }
    //std::cout<<"nodes_order_size="<<nodes_order.size()<<"\n";
    /*std::cout<<"contents of nodes_order are";
    for (std::vector<uint>::iterator it = nodes_order.begin(); it != nodes_order.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';*/

}

void Ethereum::perc_links(double perc_new){
    uint num_to_remove= perc_new*(M-num_exchange_links);
    //std::cout<<"removed="<<num_to_remove<<" for"<<perc_new<<" M="<<M<<"\n";
    for (uint uu=0; uu<num_to_remove; ++uu){
        vector<uint> link_to_remove=all_links[uu];
        link_alive[link_to_remove]=false;

    }
    //std::cout<<"nodes_order_size="<<nodes_order.size()<<"\n";
    /*std::cout<<"contents of nodes_order are";
    for (std::vector<uint>::iterator it = nodes_order.begin(); it != nodes_order.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';*/

}

void Ethereum::perc_links_special(double perc_share_share, double perc_noshare_share,double perc_share_noshare){
    for(auto link : all_links_share_share){
        link_alive[link]=true;
    }
    for(auto link : all_links_noshare_share){
        link_alive[link]=true;
    }
    for(auto link : all_links_share_noshare){
        link_alive[link]=true;
    }
    uint num_to_remove_share_share= perc_share_share*all_links_share_share.size();
    std::cout<<"remove="<<num_to_remove_share_share<<" share-share="<<all_links_share_share.size()<<"\n";
    for (uint uu=0; uu<num_to_remove_share_share; ++uu){
        vector<uint> link_to_remove=all_links_share_share[uu];
        link_alive[link_to_remove]=false;
    }
    uint num_to_remove_noshare_share= perc_noshare_share*all_links_noshare_share.size();
    for (uint uu=0; uu<num_to_remove_noshare_share; ++uu){
        vector<uint> link_to_remove=all_links_noshare_share[uu];
        link_alive[link_to_remove]=false;
    }
    uint num_to_remove_share_noshare= perc_share_noshare*all_links_share_noshare.size();
    for (uint uu=0; uu<num_to_remove_share_noshare; ++uu){
        vector<uint> link_to_remove=all_links_share_noshare[uu];
        link_alive[link_to_remove]=false;
    }
    //std::cout<<"nodes_order_size="<<nodes_order.size()<<"\n";
    /*std::cout<<"contents of nodes_order are";
    for (std::vector<uint>::iterator it = nodes_order.begin(); it != nodes_order.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';*/

}


void Ethereum::writeDist_directed(uint runs, double perc_share_share, double perc_noshare_share,double perc_share_noshare, string text){
        std::ofstream outfile3;
        std::ostringstream oss3;
    
        oss3<<"Exchange_dist_"<<N<<"_share_share_perc"<<(int) 10*perc_share_share<<"_"<<"_noshare_share_perc"<<(int) 10*perc_noshare_share<<"_share_noshare_perc"<<(int) 10*perc_share_noshare<<"_"<<text<<"_";
        oss3<<"_"<<runs;

        outfile3.open((oss3.str()+".txt").c_str());
        //outfile2<< "num_links_inter" <<"    "<<"num_links_intra"<<"    "<<"gcc\n";

        for (uint mm=0; mm<distance_exchange.size(); mm++){
            outfile3<<mm<<" "<<distance_exchange[mm]<<"\n";}
        outfile3.close();
}

void Ethereum::writeDist(uint runs, double perc_val, string text){
        std::ofstream outfile3;
        std::ostringstream oss3;
    
        oss3<<"Exchange_dist_"<<N<<"_perc"<<(int) 10*perc_val<<"_"<<text<<"_";
        oss3<<"_"<<runs;

        outfile3.open((oss3.str()+".txt").c_str());
        //outfile2<< "num_links_inter" <<"    "<<"num_links_intra"<<"    "<<"gcc\n";

        for (uint mm=0; mm<distance_exchange.size(); mm++){
            outfile3<<mm<<" "<<distance_exchange[mm]<<"\n";}
        outfile3.close();
}

void Ethereum::writeS(uint runs, string text)
{
    std::ofstream outfile2;
    std::ostringstream oss2;
    
    oss2<<"S_vals_ethereum_100comp_"<<N<<"_"<<text<<"_";
    oss2<<"_"<<runs;

    outfile2.open((oss2.str()+".txt").c_str());
    //outfile2<< "num_links_inter" <<"    "<<"num_links_intra"<<"    "<<"gcc\n";

    for (uint mm=0; mm<numlinks_S_history.size(); mm++){
        for (uint nn=0; nn<numlinks_S_history[mm].size(); nn++){
            std::cout<<numlinks_S_history[mm][nn]<<" ";
            outfile2<< numlinks_S_history[mm][nn]<<"    ";}
        std::cout<<"\n";
        outfile2<<"\n";}
    outfile2.close();
}

vector<int> Ethereum::get_path_lengths(uint node){
       uint  distval, j;
       vector<int> parent(N,-1); //save parent to get to shortest path
       pathlength=vector<int>(N, 999999999);
       
       //Step1 find shortest path from an exchange node
       std::fill(parent.begin(), parent.end(), -1); //stores the nodes parent on the shortest path
       std::fill(pathlength.begin(), pathlength.end(),999999999); //stores the current path length to the node
       //condition if we have found the shortet path yet (we do not need to find all nodes through bfs)
       //std::cout<<"known="<<known_nodes.size()<<"\n";
       std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
       while(!bfs_queue.empty()){ bfs_queue.pop();}
        bfs_queue.push(node);
        parent[node]=-1;
        pathlength[node]=0;
       while(!bfs_queue.empty()){
           j=bfs_queue.front();

           bfs_queue.pop();
           for(auto neighbor : adjacency_list[j]){
               if( bfs_visited[neighbor] == BFS::White ){
                   distval=(pathlength[j]+1);
                   if(distval<pathlength[neighbor]){
                       pathlength[neighbor]=distval;
                       parent[neighbor]=j; //The parent of this node is j if we are seeing it for the first time (BFS::White)
                   }
                   bfs_queue.push(neighbor); // <-- only path that increases queue
                   bfs_visited[neighbor] = BFS::Gray;
               

                }
            }
               
        
           //std::cout<<"conditions="<<node_spt_reached<<" "<<bfs_queue.size()<<"\n";
           bfs_visited[j]=BFS::Black;

       }
    return (pathlength);
}

void Ethereum::reach_node(uint node){
    /* Algorithm to determine how many interrogations it would actually take an
    investigator to reach a given node assuming the investigator has access to
    the exchange nodes.
    (1) From full network (ignoring the link_alive parameter), find the k shortest paths to the
    node from each exchange node.
    (2) Beginning with the shortest path:
        (a)  Check the first link from the exchange to the path, it of course exists.
        (b) Check the next link from the discovered node
            and increment total interrogations by one, if the link exists, continue.
        (c) If following all links on the path is possible (incrementing total interrogations
            each time), then the desired node is found, and that is our shortest path.
        (d) If a link along the path does not exist, set that link_alive_reach, to false.
        (e) add the discovered nodes to the (modified) exchange node list
            and find new shortest path.
        (f) Repeat steps (a) through (e) until the desired node is found.
    */
    path_overall=999999999;
    uint  distval, j;
    bool node_reached_overall=false, node_spt_reached=false, known_unreachable=false;
    map< vector< uint >, bool > link_known_alive(link_alive);
    for(auto ll : link_known_alive){
        link_known_alive[ll.first]=true; //at first we believe all links to be alive.
    }
    
    vector<int> parent(N,-1); //save parent to get to shortest path
    pathlength=vector<int>(N, 999999999);
    vector<uint> known_nodes(exchange_nodes);
    /*
    BFS
    distance_so_far=0
    while node_found==false:
        distance_so_far+=1
        for each node that is known:
            traverse all its links
            if node is not in known_nodes, add node to known_nodes
            if node is found
    */
    
    while(node_reached_overall==false and known_unreachable==false){
    
        //Step1 find shortest path from an exchange node
        std::fill(parent.begin(), parent.end(), -1); //stores the nodes parent on the shortest path
        std::fill(pathlength.begin(), pathlength.end(),999999999); //stores the current path length to the node
        //condition if we have found the shortet path yet (we do not need to find all nodes through bfs)
        //std::cout<<"known="<<known_nodes.size()<<"\n";
        std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
        node_spt_reached=false;
        while(!bfs_queue.empty()){ bfs_queue.pop();}

        for (auto i: known_nodes){
                bfs_queue.push(i);
    
            if(bfs_visited[i] == BFS::White and alive[i]){
                bfs_queue.push(i);
                bfs_visited[i]=BFS::Gray;
                pathlength[i]=0; //pathlength of an exchange to an exchange is 0
                if(std::find(exchange_nodes.begin(), exchange_nodes.end(), i) != exchange_nodes.end()){
                    parent[i]=-1;
                    }
                }
        }
        while(node_spt_reached==false and (!bfs_queue.empty())){
            j=bfs_queue.front();
            if(j==node){ //this should never really happen, but we can leave it here anyway.....
                node_spt_reached=true;
                std::cout<<"Warning, node reached in queue\n";
            }
            bfs_queue.pop();
            for(auto neighbor : adjacency_list[j]){
                
                if( bfs_visited[neighbor] == BFS::White and alive[neighbor]){
                    if(j<neighbor){
                        link1={j,neighbor};
                    }
                    
                    if(neighbor<j){
                        link1={neighbor,j};
                    }
                    
                    
                    if (link_known_alive[link1]){
                        distval=(pathlength[j]+1);
                        if(distval<pathlength[neighbor]){
                            pathlength[neighbor]=distval;
                            parent[neighbor]=j; //The parent of this node is j if we are seeing it for the first time (BFS::White)


                        }
                        bfs_queue.push(neighbor); // <-- only path that increases queue
                        bfs_visited[neighbor] = BFS::Gray;
                    
                    
                        if (neighbor==node){
                            node_spt_reached=true;
                            //std::cout<<"parent="<<parent[neighbor]<<" "<<i<<"\n";
                            //std::cout<<"parent="<<parent[neighbor]<<"\n";

                        }
                    }
                }
            }
            //std::cout<<"conditions="<<node_spt_reached<<" "<<bfs_queue.size()<<"\n";
            bfs_visited[j]=BFS::Black;

        }
        
    
        //std::cout<<"believed_pathlength="<<pathlength[node]<<"\n";
        //std::cout<<node<<" parent="<<parent[node]<<"\n";

        /*uint node_found=node;
        std::cout<<node<<"'s path=";
        while(std::find(exchange_nodes.begin(), exchange_nodes.end(), node_found) == exchange_nodes.end()){
            node_found=parent[node_found];
            std::cout<<node_found<<" ";

        }
        std::cout<<"\n";*/
        
        //Now attempt to traverse the shortest path
        bool path_checked=false;
        bool path_works=true;
        uint current_node=node;
        vector<uint> nodes_discovered_now;
        vector<uint>  link_to_kill={0,0};
        int k;
        while(path_checked==false){
            k=parent[current_node]; //need here a uint, but later if the node exists we must convert to uint
            //std::cout<<current_node<<" parent="<<k<<"\n";
            if(k==-1 or known_nodes.size()-exchange_nodes.size()>25){ //I cannot find any path to this node now.
                path_checked=true;
                known_unreachable=true;
                path_works=false;
                break;
            }
            j=k; //here to convert to uint
            //std::cout<<"node on="<<j<<"\n";
            if(j<current_node){
                link1={j,current_node};
            }

            if(current_node<j){
                link1={current_node,j};
            }
            if(link_alive[link1]){ //if the link is actually alive
                nodes_discovered_now.push_back(current_node);
                //std::cout<<"added to discovered="<<current_node<<"\n";
                }

            if(link_alive[link1]==false){ //here we check the actual link_alive
                //std::cout<<"failed link="<<link1[0]<<" "<<link1[1]<<"\n";
                link_to_kill=link1; //This is the link we will discover is not known
                path_works=false; //we now will know that the path fails since one link is false
                nodes_discovered_now.clear(); //we did not discover any nodes since we are going backwards, we will discover all prior nodes to the failed link though
            }

            if(std::find(known_nodes.begin(), known_nodes.end(), j) != known_nodes.end()){
                path_checked=true;
            }
            current_node=j;
        }
        if(path_works==true){ //if my last path had no link_alive==false then I have discovered the node
            node_reached_overall=true;
        }
        if(path_works==false and known_unreachable==false){
                link_known_alive[link1]=false;
        }
        //std::cout<<"discovered= ";
        for(auto d : nodes_discovered_now){
            if(std::find(known_nodes.begin(), known_nodes.end(), d) == known_nodes.end()){
                //std::cout<<d<<" ";
                known_nodes.push_back(d);
                for (auto neighbor : adjacency_list[d] ){
                        if(d<neighbor){
                            link1={d,neighbor};
                        }
                    
                        if(neighbor<d){
                            link1={neighbor,d};
                        }
                        link_known_alive[link1]=link_alive[link1];
                }
            }
        }
        //std::cout<<"\n";

        if(path_works==true){
            path_overall=(int) known_nodes.size()- (int) exchange_nodes.size();
        }
        if(path_works==false and known_unreachable==true){ //still save this number of interrogations, but also note that path was not found by making it negative
            path_overall=-((int) known_nodes.size()-(int) exchange_nodes.size());
        }

    }
    //std::cout<<"path_overall="<<path_overall<<"\n";
}












void Ethereum::BFS1()
{

    std::fill(components.begin(),components.end(),-1);
    std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
    while(!bfs_queue.empty()){ bfs_queue.pop();}
    component_size.clear();
    uint j,distval=0,counter_exchange=0,compsize=0;
    uint num_exchange_nodes=exchange_nodes.size();
    bool is_exchange;
    for(auto i : bfs_order_exchanges){
        is_exchange=false;
        if(counter_exchange<num_exchange_nodes) is_exchange=true;
        compsize=0;
        counter_exchange++;
        if(is_exchange){
            std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
            component_size.clear();
            while(!bfs_queue.empty()){ bfs_queue.pop();}
            std::fill(components.begin(),components.end(),-1);
            if(bfs_visited[i] == BFS::White and alive[i]){
                bfs_queue.push(i);
                bfs_visited[i]=BFS::Gray;
                components[i]=i;
                compsize++;
                distance_exchange[i]=0;
                while(!bfs_queue.empty()){
                    j=bfs_queue.front();
                    bfs_queue.pop();
                    for(auto neighbor : adjacency_list[j]){
                        if( bfs_visited[neighbor] == BFS::White and alive[neighbor]){
                            if(j<neighbor){
                                link1={j,neighbor};
                            }
                            if(neighbor<j){
                                link1={neighbor,j};
                            }
                            if (link_alive[link1]){
                                distval=(distance_exchange[j]+1);
                                if(distval<distance_exchange[neighbor]){
                                    distance_exchange[neighbor]=distval;
                                }
                                bfs_queue.push(neighbor); // <-- only path that increases queue
                                bfs_visited[neighbor] = BFS::Gray;
                                components[neighbor] = i;
                                compsize++;
                            }
                        }
                    }
                    bfs_visited[j]=BFS::Black;

                }
            }
            component_size[i] = compsize;
        }
        
        else{
            if(bfs_visited[i] == BFS::White and alive[i]){
                bfs_queue.push(i);
                bfs_visited[i]=BFS::Gray;
                components[i]=i;
                compsize++;
                while(!bfs_queue.empty()){
                    j=bfs_queue.front();
                    bfs_queue.pop();
                    for(auto neighbor : adjacency_list[j]){
                        if( bfs_visited[neighbor] == BFS::White and alive[neighbor]){
                            if(j<neighbor){
                                link1={j,neighbor};
                            }
                            if(neighbor<j){
                                link1={neighbor,j};
                            }
                            if (link_alive[link1]){
                                bfs_queue.push(neighbor); // <-- only path that increases queue
                                bfs_visited[neighbor] = BFS::Gray;
                                components[neighbor] = i;
                                compsize++;
                            }
                        }
                    }
                    bfs_visited[j]=BFS::Black;

                }
            }
            component_size[i] = compsize;
        }
    }
}


void Ethereum::BFS1_directed()
{

    std::fill(components.begin(),components.end(),-1);
    std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
    while(!bfs_queue.empty()){ bfs_queue.pop();}
    component_size.clear();
    uint j,distval=0,counter_exchange=0,compsize=0;
    uint num_exchange_nodes=exchange_nodes.size();
    bool is_exchange;
    for(auto i : bfs_order_exchanges){
        is_exchange=false;
        if(counter_exchange<num_exchange_nodes) is_exchange=true;
        compsize=0;
        counter_exchange++;
        if(is_exchange){
            std::fill(bfs_visited.begin(), bfs_visited.end(),BFS::White);
            component_size.clear();
            while(!bfs_queue.empty()){ bfs_queue.pop();}
            std::fill(components.begin(),components.end(),-1);
            if(bfs_visited[i] == BFS::White and alive[i]){
                bfs_queue.push(i);
                bfs_visited[i]=BFS::Gray;
                components[i]=i;
                compsize++;
                distance_exchange[i]=0;
                while(!bfs_queue.empty()){
                    j=bfs_queue.front();
                    bfs_queue.pop();
                    for(auto neighbor : adjacency_list[j]){
                        if( bfs_visited[neighbor] == BFS::White and alive[neighbor]){
                            link1={neighbor,j};
                            if (link_alive[link1]){
                                distval=(distance_exchange[j]+1);
                                if(distval<distance_exchange[neighbor]){
                                    distance_exchange[neighbor]=distval;
                                }
                                bfs_queue.push(neighbor); // <-- only path that increases queue
                                bfs_visited[neighbor] = BFS::Gray;
                                components[neighbor] = i;
                                compsize++;
                            }
                        }
                    }
                    bfs_visited[j]=BFS::Black;

                }
            }
            component_size[i] = compsize;
        }
        
        else{
            if(bfs_visited[i] == BFS::White and alive[i]){
                bfs_queue.push(i);
                bfs_visited[i]=BFS::Gray;
                components[i]=i;
                compsize++;
                while(!bfs_queue.empty()){
                    j=bfs_queue.front();
                    bfs_queue.pop();
                    for(auto neighbor : adjacency_list[j]){
                        if( bfs_visited[neighbor] == BFS::White and alive[neighbor]){
                            link1={neighbor,j};
                            if (link_alive[link1]){
                                bfs_queue.push(neighbor); // <-- only path that increases queue
                                bfs_visited[neighbor] = BFS::Gray;
                                components[neighbor] = i;
                                compsize++;
                            }
                        }
                    }
                    bfs_visited[j]=BFS::Black;

                }
            }
            component_size[i] = compsize;
        }
    }
}


