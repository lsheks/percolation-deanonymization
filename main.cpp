//
//  main.cpp
//  ethereum
//
//  Created by Louis Shekhtman on 9/13/18.
//  Copyright © 2018 Louis Shekhtman. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "ethereum.cpp"
#include <math.h>
#include <ctime>
using std::vector;


int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;

    uint runstart= argc > 1 ? atoi(argv[1]) : 0;
    uint runs= argc > 2 ? atoi(argv[2]) : 1;

    uint steps=5;
    //double dom_val=(double) 0.7;
    //std::srand(time(0));

    vector<double> perc_to_sample1;
    //perc_to_sample1.push_back(0.);
    perc_to_sample1.push_back(0.8);
    //perc_to_sample1.push_back(0.9);
    //perc_to_sample1.push_back(0.98);
    //perc_to_sample1.push_back(0.99);
    //perc_to_sample1.push_back(1);

    /*for(uint i=0;i<steps+1;++i){
        double val=(double) (i)/(steps);
        perc_to_sample1.push_back(val);
    }*/
    
    vector<double> perc_to_sample2;
    //perc_to_sample2.push_back(0.);
    //perc_to_sample2.push_back(0.8);
    perc_to_sample2.push_back(0.9);
    perc_to_sample2.push_back(0.99);
    perc_to_sample2.push_back(1);

    uint steps2=5;
    vector<double> perc_to_sample3;
    perc_to_sample3.push_back(0.);
    perc_to_sample3.push_back(0.8);
    perc_to_sample3.push_back(0.9);
    perc_to_sample3.push_back(0.98);
    perc_to_sample3.push_back(0.99);
    perc_to_sample3.push_back(1);

    /*for(uint i=0;i<steps2+1;++i){
        double val=(double) (i)/(steps2);
        perc_to_sample3.push_back(val);
    }*/
    
    vector<uint> nodes_to_check;
    nodes_to_check.push_back(7104);
    nodes_to_check.push_back(2405);
    nodes_to_check.push_back(1247);
    //Nodes not sharers
    nodes_to_check.push_back(1077);
    nodes_to_check.push_back(9647);
    nodes_to_check.push_back(301);
    nodes_to_check.push_back(3238);
    nodes_to_check.push_back(9889);
    nodes_to_check.push_back(1224);
    nodes_to_check.push_back(3755);

    #pragma omp parallel for num_threads(5)
    for (uint node_choice=0; node_choice<10; node_choice++){
        uint ll=nodes_to_check[node_choice];//rand() % 10407;
        for( uint run=runstart; run<runs; run++){
            std::srand(time(0)+run);
            Ethereum Ethereum1(run);
            //std::cout<<k_vals.size()<<" "<<m_mod.size()<<"\n";
            Ethereum1.read_nodes_dacunha_directed();
            uint exchange_node= nodes_to_check[node_choice]; //rand() % 10407;
            vector<uint> exchange_nodes={exchange_node};
            Ethereum1.link_order_matjaz_dacunha(exchange_nodes);
            //Ethereum1.order_by_deg(); //For ordering the nodes by degree


            for(auto perc_val1 : perc_to_sample1){
                for(auto perc_val2 : perc_to_sample2){
                    for(auto perc_val3 : perc_to_sample3){

                        //Ethereum1.perc_nodes(perc_val);
                        Ethereum1.perc_links_special(perc_val1, perc_val2, perc_val3);
                        std::cout<<"perc1="<<perc_val1<<", perc2="<<perc_val2<<", perc3="<<perc_val3<<", run="<<run<<"\n";
                        Ethereum1.update_S_directed(1-perc_val1); //take 1 minus to get survival prob
                        if( ((int) (100*perc_val1) % 1)==0){
                            string add_text="da_cunha_other_search_biased_";
                            add_text+=std::to_string(exchange_node);
                            add_text+="_";
                            Ethereum1.writeDist_directed(run, 1-perc_val1,1-perc_val2,1-perc_val3, add_text);
                        }
                        string add_text1="da_cunha_other_search_biased_";
                        add_text1+=std::to_string(exchange_node);
                        add_text1+="_";
                        Ethereum1.writeS(run, add_text1); //for treelike
                    }
                }

            }

        }
    }
    return 0;
}


