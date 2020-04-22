#include <set>
#include <iostream>
#include <networkit/graph/Graph.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <bits/stdc++.h>
#include<map>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <armadillo>
#include <matrixPreparation.h>
namespace po = boost::program_options;
using namespace std::chrono;
using namespace std;
using namespace arma;


void readVector(string filename, vector<string>& geneNode){
    std::string line;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        getline (myfile,line);

        boost::algorithm::split(geneNode, line, boost::is_any_of(","));
    }
}

void readVector(string filename, vector<NetworKit::node>& goSelected){
    std::string line;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        getline (myfile,line);
        vector<string> res;
        boost::algorithm::split(res, line, boost::is_any_of(","));

        for(string s : res){
            NetworKit::node value;
            std::istringstream iss(s);
            iss >> value;
            goSelected.push_back(value);
        }

    }
}



mat computeR(mat& A, mat& B, mat& annotation, double alpha){
    return alpha* (A*B) + (1.-alpha)*annotation;
}

void runRW(mat& hybrid, mat& tpMatrix, mat& annotation,mat& predictionMatrix,int iter_P,int iter_F, double alpha){

    mat RpOld = annotation;
    mat  RfOld = annotation;
    map<int,mat> iter2mat;
    int maxIter = max(iter_P,iter_F);
    for(int iter=1; iter<=maxIter;iter++){
        if(iter<=iter_P && iter <= iter_F){
            mat  Rf = computeR(RfOld,tpMatrix,annotation,alpha);
            mat  Rp = computeR(hybrid,RpOld,annotation,alpha);
            //  mat  Rf = computeR(RfOld,tpMatrix,RfOld,alpha,"RF",iter);
            //  mat  Rp = computeR(hybrid,RpOld,RpOld,alpha,"RP",iter);
            mat newIt = (Rf+Rp)/2.;
            iter2mat.insert(pair<int,mat>(iter,newIt));
            RpOld = newIt;//Rp;
            RfOld = newIt;//Rf;

        }else if(iter<=iter_P){
            mat   Rp = computeR(hybrid,RpOld,annotation,alpha);

            mat newIt = (RfOld+Rp)/2.;
            iter2mat.insert(pair<int,mat>(iter,newIt));
            RpOld = Rp;
        }else if(iter<=iter_F){
            mat Rf = computeR(RfOld,tpMatrix,annotation,alpha);
            mat newIt = (Rf+RpOld)/2.;
            iter2mat.insert(pair<int,mat>(iter,newIt));
            RfOld = Rf;
        }
    }
    predictionMatrix = iter2mat.at(maxIter);
}



int newGOA(int ac, char* av[])
{


    int threads;
    string inFolder;
    string outfile;
    string obofile;
    int iter_P;
    int iter_F;
    double alpha;


    po::options_description desc(" options");
    desc.add_options()
            ("help,h", "produce help message")
            ("num_threads,T", po::value<int>(&threads) ->default_value(1), "set number of threads to parallize")
            ("inputFolder I", po::value<string>(&inFolder)->default_value("./results"), "Set an input folder. This folder was created with the matrixPreparation function.")
            ("goFile,G",po::value<string>(&obofile),"Set the path of gene ontology file [OBO format]")
            ("hybrid,H","Computing newGOA with the hybrid matrix (ppi + semanticPpi)")
            ("outFileName,O",po::value<string>(&outfile)->default_value("predictedAnnotation.txt"), "output file name to provide the results")
            ("iter_P,P",po::value<int>(&iter_P)->default_value(5),"Maximal number of iteration to the ppi network")
            ("iter_F,F",po::value<int>(&iter_F)->default_value(5),"Maximal number of iteration to the go network")
            ("alpha,A", po::value<double>(&alpha)->default_value(0.1), "set the probability for a walker staying at the starting point")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    po::store(po::command_line_parser(ac, av).
              options(desc).run(), vm);
    //  cout<<ac<<endl;
    if (vm.count("help") || ac == 2)
    {
        std::cerr << desc << '\n';
        return 1;
    }
    vector<string> err;



    bool flag = false;
    if(!vm.count("goFile")){
        err.push_back("[ERROR] Please set the gene ontology file [OBO format] \n");
        flag=true;
    }


    if(flag){
        for(string s : err){
            cerr<<s;
        }
        cerr<< desc<<endl;
        return 1;
    }

    if (!boost::filesystem::exists(inFolder)){
        err.push_back("[ERROR] Please use a correct path including the matrixPreparation result folder.");
        flag=true;
    }


    if (!boost::filesystem::exists(obofile)){
        err.push_back("[ERROR] Please use a correct path including the obo file.");
        flag=true;
    }



    if(flag){
        for(string s : err){
            cerr<<s;
        }
        cerr<< desc<<endl;
        return 1;
    }


    NetworKit::Graph go(0,true,true); // the first true is used to consider the graph as weigthed and the second true to consider the graph directed
    map<NetworKit::node,GOTerm> idx2goterm;
    map<string,NetworKit::node> gotermsID2idx;
    map<string,Root> namespace2root;

    auto      start = high_resolution_clock::now();
    reader(obofile,idx2goterm, go, gotermsID2idx,namespace2root,threads);
    auto  stop = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(stop - start);
    cout <<"To read the Gene Ontology took " <<duration2.count()<<" milliseconds." << endl;
    start = high_resolution_clock::now();
    vector<string> geneNode;
    readVector(inFolder+"/genes.txt", geneNode);
    vector<NetworKit::node> goSelected;
    readVector(inFolder+"/GOterms.txt", goSelected);
    stop = high_resolution_clock::now();
    duration2 = duration_cast<milliseconds>(stop - start);
    cout <<"To read the gene and GO term node files took " <<duration2.count()<<" milliseconds." << endl;
    start = high_resolution_clock::now();
    mat proteinPPI;
    proteinPPI.load(inFolder+"/nPpi.bin");
    mat semanticPPI;
    semanticPPI.load(inFolder+"/nSSPpi.bin");
    mat hybridPPI;
    hybridPPI.load(inFolder+"/hybrid.bin");
    cout<<"start"<<endl;
    mat TP;
    TP.load(inFolder+"/TPMatrix.bin");
    //sp_mat sparseTP(TP);



    mat annotation;
    annotation.load(inFolder+"/annotationMatrix.bin");
    stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout <<"To read the matrix files (ppis and sparse TP) took " <<duration.count()<<" seconds." << endl;

    // READ GO file

    start = high_resolution_clock::now();
    mat predictionMatrix;
    if(vm.count("hybrid")){
        runRW(hybridPPI,TP,annotation,predictionMatrix,iter_P,iter_F,alpha);
    }else{
        runRW(proteinPPI,TP,annotation,predictionMatrix,iter_P,iter_F,alpha);
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    cout <<"To run the random walker took " <<duration.count()<<" seconds." << endl;

    stringstream ss;
    ss<<"gene,GOterm,OriginalAnnotation,PredictionScore\n";
    for(size_t i=0; i<predictionMatrix.n_rows;i++){

        for(size_t j = 0; j< predictionMatrix.n_cols; j++){
            if(predictionMatrix(i,j) > 0.0){

                double pv = predictionMatrix(i,j);
                double av = annotation(i,j);
                ss<<geneNode.at(i)<<","<<idx2goterm.at(goSelected.at(j)).id<<","<<av<<","<<pv<<"\n";
            }
        }
    }
    ofstream res;
    res.open(inFolder+"/"+outfile);
    res<< ss.str();
    res.close();

    return 0;
}


