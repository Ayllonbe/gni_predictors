#pragma once
#ifndef MATRIXPREPARATION_H
#define MATRIXPREPARATION_H

#include <set>
#include <iostream>
#include <networkit/graph/Graph.hpp>
#include <bits/stdc++.h>
#include<map>
#include <robin_hood.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/minmax.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <armadillo>

const NetworKit::edgeweight ew_default = 1.0;
// GO Term class
// That include the information of a GO term have (only the interesting ones)
// id - the GOID (format GO:XXXXX)
// name - the GO name
// top - The GO namespace (biological_process, molecular_function or cellular_component)
// index - the node index in the created graph
// private ic - the information content
class GOTerm{
public:
    std::string id;
    std::string name;
    std::string top;
    NetworKit::node index; // this index corresponds to graph node index.
    GOTerm(){
    }
    GOTerm(std::string& id, std::string& name, std::string& top,NetworKit::node& index){
        this->id=id;
        this->name=name;
        this->top=top;
        this->index=index;
    }
    GOTerm(std::vector<std::string>& termInfo,NetworKit::node& index){
        this->id=termInfo.at(0);
        this->name=termInfo.at(1);
        this->top=termInfo.at(2);
        this->index=index;
    }
    void set_IC(double& val){
        this->ic=val;
    }
    double get_IC(){
        return this->ic;
    }
    void set_child(NetworKit::node c){
        this->childs.push_back(c);
    }
    NetworKit::node get_child(int pos){
        return this->childs.at(pos);
    }
    std::vector<NetworKit::node>& get_childrens(){
        return this->childs;
    }
    void set_parent(NetworKit::node c){
        this->parents.push_back(c);
    }
    NetworKit::node& get_parent(int pos){
        return this->parents.at(pos);
    }
    std::vector<NetworKit::node>& get_parents(){
        return this->parents;
    }

private:
    double ic;
    std::vector<NetworKit::node> childs;
    std::vector<NetworKit::node> parents;
    //     vector<NetworKit::node> descendants;
    //      vector<NetworKit::node> ancestors;

};


robin_hood::unordered_set<NetworKit::node> getdescNodes(NetworKit::Graph& go, NetworKit::node& n);


std::vector<GOTerm> getdescendants(NetworKit::Graph& go, NetworKit::node& n, robin_hood::unordered_map<NetworKit::node,GOTerm>& idx2goterm);

robin_hood::unordered_set<NetworKit::node> getancNodes( NetworKit::node& n ,robin_hood::unordered_map<NetworKit::node,GOTerm>& idx2goterm);


std::vector<GOTerm> getancestors(NetworKit::node& n, robin_hood::unordered_map<NetworKit::node,GOTerm>& idx2goterm);


class Root{
public:
    std::string id;
    long ndescendants;
    Root(std::string& id, long& ndesc){
        this->id=id;
        this->ndescendants=ndesc;
    }

};



// Reading obo file
// This function allow to read a go file and fulfill three variables:
//      idx2goterm - a map with the node index as key and the GOTerm class as
//      go - Graph of Gene Ontology
//      gotermsIDidx - a map with the GOID as key and the node index as value
void reader(std::string go_file, robin_hood::unordered_map<NetworKit::node,GOTerm>& idx2goterm,NetworKit::Graph& go, robin_hood::unordered_map<std::string,NetworKit::node>& gotermsID2idx, robin_hood::unordered_map<std::string,Root>& namespace2root, int& threads);

class Annotation{
public:
    Annotation(){

    }
    Annotation(robin_hood::unordered_map<NetworKit::node, std::vector<std::string>>& go2gvec, robin_hood::unordered_map<std::string, robin_hood::unordered_set<NetworKit::node>>& gene2govec){
        this->go2gvec = go2gvec;
        this->gene2govec =gene2govec;
    }
    void setAnnotation(NetworKit::node go, robin_hood::unordered_set<std::string>& genes){

        std::vector<std::string> vec(genes.begin(),genes.end());
        this->go2gvec.insert({go,vec});

        for(std::string g : genes){
            this->go2gvec.at(go).push_back(g);
            if(this->gene2govec.find(g)==this->gene2govec.end()){
                robin_hood::unordered_set<NetworKit::node> vec;
                this->gene2govec.insert({g,vec});
            }
            this->gene2govec.at(g).insert(go);
        }
    }
    std::vector<std::string> getGenes(NetworKit::node& goTerm){
        if(this->go2gvec.find(goTerm) != this->go2gvec.end()){
            return this->go2gvec.at(goTerm);
        }else{
            return std::vector<std::string>();
        }

    }
    robin_hood::unordered_set<NetworKit::node>& getGOTerms(std::string& gene){
        return this->gene2govec.at(gene);
    }
    int findGene(std::string& gene){
        return this->gene2govec.find(gene)!=gene2govec.end();
    }
    int findGO(NetworKit::node& go){
        return this->go2gvec.find(go)!=go2gvec.end();
    }
private:
    robin_hood::unordered_map<NetworKit::node, std::vector<std::string>> go2gvec;
    robin_hood::unordered_map<std::string, robin_hood::unordered_set<NetworKit::node>> gene2govec;
};



double icSum(std::vector<NetworKit::node>& childrens, robin_hood::unordered_map<NetworKit::node, GOTerm>& idx2goTerm);

int matrixPreparation(int argc, char * argv[]);

#endif // MATRIXPREPARATION_H
