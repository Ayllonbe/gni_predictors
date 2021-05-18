#include <matrixPreparation.h>
namespace po = boost::program_options;
using namespace std::chrono;

robin_hood::unordered_set<NetworKit::node> getdescNodes(NetworKit::Graph &go,
                                                        NetworKit::node &n) {
    robin_hood::unordered_set<NetworKit::node> descendants;
    NetworKit::Graph::NeighborRange<true> childIT = go.inNeighborRange(n);

    for (NetworKit::Graph::NeighborIterator it = childIT.begin();
         it != childIT.end(); it++) {
        NetworKit::node nd = it.operator*();
        descendants.insert(nd);
        robin_hood::unordered_set<NetworKit::node> vd = ::getdescNodes(go, nd);
        descendants.insert(vd.begin(), vd.end());
    }
    return descendants;
};

std::vector<GOTerm>
getdescendants(NetworKit::Graph &go, NetworKit::node &n,
               robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goterm) {
    std::vector<GOTerm> descendants;
    robin_hood::unordered_set<NetworKit::node> setNodes = ::getdescNodes(go, n);
    for (NetworKit::node n : setNodes) {
        descendants.push_back(idx2goterm.at(n));
    }
    return descendants;
};

robin_hood::unordered_set<NetworKit::node>
getancNodes(NetworKit::node &n,
            robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goterm) {
    robin_hood::unordered_set<NetworKit::node> ancestors;

    for (NetworKit::node nd : idx2goterm.at(n).get_parents()) {

        ancestors.insert(nd);
        robin_hood::unordered_set<NetworKit::node> vd =
                ::getancNodes(nd, idx2goterm);
        ancestors.insert(vd.begin(), vd.end());
    }
    return ancestors;
};

std::vector<GOTerm>
getancestors(NetworKit::node &n,
             robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goterm) {
    std::vector<GOTerm> descendants;
    for (NetworKit::node n : ::getancNodes(n, idx2goterm)) {
        descendants.push_back(idx2goterm.at(n));
    }
    return descendants;
};

// Reading obo file
// This function allow to read a go file and fulfill three variables:
//      idx2goterm - a map with the node index as key and the GOTerm class as
//      go - Graph of Gene Ontology
//      gotermsIDidx - a map with the GOID as key and the node index as value
void reader(
        std::string go_file,
        robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goterm,
        NetworKit::Graph &go,
        robin_hood::unordered_map<std::string, NetworKit::node> &gotermsID2idx,
        robin_hood::unordered_map<std::string, Root> &namespace2root,
        int &threads) {
    std::string line;
    std::ifstream myfile(go_file);
    std::string delimiter = ":";
    std::string token = "";
    bool isTerm = false;
    std::vector<std::string> termInfo;
    bool obsolete = false;
    robin_hood::unordered_map<NetworKit::node, std::vector<std::string>>
            node2parent;
    std::vector<std::string> parentvector;
    std::vector<NetworKit::node> rootID;

    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            if (line == "") {
                token = line;
                if (isTerm == true && obsolete != 1) {

                    NetworKit::node n = go.addNode();
                    GOTerm term(termInfo, n);
                    idx2goterm.insert({n, term});

                    gotermsID2idx.insert({termInfo.at(0), n});

                    if (parentvector.size() == 0) {
                        rootID.push_back(n);
                    }

                    node2parent.insert({n, parentvector});
                }
                isTerm = 0;
                parentvector.clear();
                termInfo.clear();
                isTerm = false;
                obsolete = false;
            }

            if (line == "[Term]") {
                isTerm = 1;
                token = line;
            }

            if (token == "[Term]") {

                std::string key = line.substr(0, line.find(delimiter));
                std::string value = line.substr(line.find(delimiter) + 2, line.size());
                // cout << key <<"-:-"<<value<<"\n";
                if (key == "id") {
                    termInfo.push_back(value);
                } else if (key == "name") {
                    termInfo.push_back(value);
                } else if (key == "is_a") {
                    std::string delimiter = " ! ";
                    std::string rel = value.substr(0, value.find(delimiter));
                    parentvector.push_back(rel);

                } else if (key == "namespace") {
                    termInfo.push_back(value);
                } else if (key == "is_obsolete") {
                    if (value == "true") {
                        obsolete = true;
                    }
                }
            }
        }
        myfile.close();

        robin_hood::unordered_map<NetworKit::node,
                std::vector<std::string>>::iterator it;
        for (it = node2parent.begin(); it != node2parent.end(); it++) {

            for (std::string v : it->second) {

                go.addEdge(it->first, gotermsID2idx.at(v), ew_default);
                idx2goterm.at(it->first).set_parent(gotermsID2idx.at(v));
                idx2goterm.at(gotermsID2idx.at(v)).set_child(it->first);
            }
        }
        for (NetworKit::node n : rootID) {
            long nd = static_cast<long>(::getdescendants(go, n, idx2goterm).size());
            GOTerm gt = idx2goterm.at(n);
            Root r(gt.id, nd);
            namespace2root.insert({gt.name, r});
        }

        robin_hood::unordered_map<NetworKit::node, GOTerm>::iterator iter =
                idx2goterm.begin();
        // #pragma omp parallel num_threads(threads)
        //     {

        for (; iter != idx2goterm.end(); iter++) {
            NetworKit::node n = iter->first;
            GOTerm &gt = iter->second;
            double ic =
                    1 - log10(static_cast<double>(::getdescNodes(go, n).size()) + 1) /
                    log10(static_cast<double>(
                              namespace2root.at(gt.top).ndescendants));
            gt.set_IC(ic);
        }
        // }

    } else {
        std::cout << "Impossible to open\n";
    }
}

arma::Mat<double> createPPIfromFile(std::string ppi_file,
                                    std::vector<std::string> &geneNodes) {
    std::ifstream myfile(ppi_file);
    std::string line;
    robin_hood::unordered_set<std::string> nodes;
    std::vector<std::string> wedges;

    if (myfile.is_open()) {
        while (getline(myfile, line)) {

            std::vector<std::string> results;

            boost::algorithm::split(results, line, boost::is_any_of("\t "));

            if (results.at(0) != "Source") {

                std::string s = results.at(0).substr(0, results.at(0).find("."));
                std::string t = results.at(1).substr(0, results.at(1).find("."));

                nodes.insert(s);
                nodes.insert(t);
                std::string l = s + " " + t + " " + results.at(2);
                wedges.push_back(l);
            }
        }
        myfile.close();

        robin_hood::unordered_map<std::string, size_t> str2id;
        size_t id = 0;
        for (std::string n : nodes) {
            geneNodes.push_back(n);
            str2id.insert({n, id});
            id++;
        }

        arma::Mat<double> ppi(geneNodes.size(), geneNodes.size(),
                              arma::fill::zeros);

        for (std::string e : wedges) {
            std::vector<std::string> results;
            boost::algorithm::split(results, e, boost::is_any_of(" "));
            double ew = ::atof(results.at(2).c_str());

            ppi.at(str2id.at(results.at(0)), str2id.at(results.at(1))) = ew;
            ppi.at(str2id.at(results.at(1)), str2id.at(results.at(0))) = ew;
        }

        return ppi;

    } else {
        std::cout << "Impossible to open\n";
    }
}

/*
 * readAnnotation function read an annotation file to create a Annotation obj.
 * This annotation obj associate a gene to the terms in the file and all their
 * ancestors. For now, take so much time, it has to be improved.
 */
void readAnnotation(
        std::string fileAnnotation, Annotation &annotation,
        robin_hood::unordered_map<std::string, NetworKit::node> &go2idx,
        robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goterm,
        std::vector<NetworKit::node> &selectedGOs) {
    std::ifstream myfile(fileAnnotation);
    std::string line;
    robin_hood::unordered_set<NetworKit::node> setGOs;

    robin_hood::unordered_map<std::string,
            robin_hood::unordered_set<NetworKit::node>>
            t2anc;
    robin_hood::unordered_map<NetworKit::node, GOTerm>::iterator gtIT;

    for (gtIT = idx2goterm.begin(); gtIT != idx2goterm.end(); gtIT++) {
        NetworKit::node n(gtIT->first);
        std::string &id = gtIT->second.id;
        t2anc.insert({id, ::getancNodes(n, idx2goterm)});
    }

    robin_hood::unordered_map<std::string, robin_hood::unordered_set<std::string>>
            go2genes;
    robin_hood::unordered_set<std::string> gotermset;
    if (myfile.is_open()) {
        while (getline(myfile, line)) {

            std::vector<std::string> results;

            boost::algorithm::split(results, line, boost::is_any_of("\t "));
            std::vector<std::string> rem;
            boost::algorithm::split(rem, results.at(0), boost::is_any_of("."));
            if (go2genes.find(results.at(1)) == go2genes.end()) {
                go2genes.insert(
                {results.at(1), robin_hood::unordered_set<std::string>{}});
                gotermset.insert(results.at(1));
            }
            go2genes.at(results.at(1)).insert(rem.at(0));
        }

        myfile.close();

        for (std::string got : gotermset) {
            if (go2idx.find(got) != go2idx.end()) {
                setGOs.insert(go2idx.at(got));

                for (NetworKit::node n : t2anc.at(got)) {
                    if (go2genes.find(idx2goterm.at(n).id) == go2genes.end()) {
                        go2genes.insert({idx2goterm.at(n).id,
                                         robin_hood::unordered_set<std::string>{}});
                        setGOs.insert(n);
                    }
                    go2genes.at(idx2goterm.at(n).id)
                            .insert(go2genes.at(got).begin(), go2genes.at(got).end());
                }
            }
        }
        robin_hood::unordered_map<
                std::string, robin_hood::unordered_set<std::string>>::iterator git;
        for (git = go2genes.begin(); git != go2genes.end(); git++) {
            std::string got(git->first);
            robin_hood::unordered_set<std::string> vec(git->second);
            if (go2idx.find(got) != go2idx.end()) {
                annotation.setAnnotation(go2idx.at(got), vec);
            }
        }

        for (NetworKit::node n : setGOs) {
            selectedGOs.push_back(n);
        }

    } else {
        std::cout << "Impossible to open\n";
    }
}

arma::Mat<double> computeGxTMatrix(std::vector<std::string> &geneNodes,
                                   std::vector<NetworKit::node> &selectedGOs,
                                   Annotation &annotation) {
    arma::Mat<double> selAssoMatrix(geneNodes.size(), selectedGOs.size(),
                                    arma::fill::zeros);
    for (size_t g = 0; g < geneNodes.size(); g++) {

        if (annotation.findGene(geneNodes.at(g))) {
            for (NetworKit::node t : annotation.getGOTerms(geneNodes.at(g))) {
                auto n = find(selectedGOs.begin(), selectedGOs.end(), t);
                auto ti = static_cast<size_t>(distance(selectedGOs.begin(), n));
                selAssoMatrix(g, ti) = 1.;
            }
        }
    }
    return selAssoMatrix;
}

double icSum(std::vector<NetworKit::node> &childrens,
             robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goTerm) {
    double sumChildsIC = 0.;

    for (NetworKit::node c : childrens) {
        sumChildsIC = sumChildsIC + idx2goTerm.at(c).get_IC();
    }
    return sumChildsIC;
}

double
transitionalProb(NetworKit::node &s, NetworKit::node &t, Annotation &annotation,
                 robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goTerm,
                 double icSum) {

    return static_cast<double>(annotation.getGenes(s).size()) /
            static_cast<double>(annotation.getGenes(t).size()) +
            idx2goTerm.at(s).get_IC() / icSum;
}

double transitionalProbSUM(
        NetworKit::node &t, Annotation &annotation,
        robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goTerm,
        double icSum) {

    double sumChildsTP = 0.;

    for (NetworKit::node c : idx2goTerm.at(t).get_childrens()) {
        sumChildsTP =
                sumChildsTP + transitionalProb(c, t, annotation, idx2goTerm, icSum);
    }

    return sumChildsTP;
}

/*
 * computeTPMatrix function computes a transational probability matrix. This
 * matrix is T x T where T includes the terms (and their ancestors) annotating
 * the genes of the chosen organism. This matrix only have a value if a term i
 * in T is associated to a term j in T.
 */
/*arma::Mat<double> computeTPMatrix_save(Annotation &annotation,
                robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goTerm,
                std::vector<NetworKit::node> &selectedGOs) {

    arma::Mat<double> TPMatrix(selectedGOs.size(), selectedGOs.size(),
                               arma::fill::zeros);
    size_t t;
    robin_hood::unordered_map<NetworKit::node, double> t2TPsum;
    robin_hood::unordered_map<NetworKit::node, double> t2ICchilds;

    for (t = 0; t < selectedGOs.size(); t++) {
        double d = icSum(idx2goTerm.at(selectedGOs.at(t)).get_childrens(), idx2goTerm);
        t2TPsum.insert(
        {selectedGOs.at(t),
         transitionalProbSUM(selectedGOs.at(t), annotation, idx2goTerm, d)});
        t2ICchilds.insert({selectedGOs.at(t), d});
    }
    for (t = 0; t < selectedGOs.size(); t++) {

        for (NetworKit::node cn :
             idx2goTerm.at(selectedGOs.at(t)).get_childrens()) {
            if (find(selectedGOs.begin(), selectedGOs.end(), cn) !=
                    selectedGOs.end()) {
                size_t pos = static_cast<size_t>(
                            distance(selectedGOs.begin(),
                                     find(selectedGOs.begin(), selectedGOs.end(), cn)));
                double tp =
                        transitionalProb(cn, selectedGOs.at(t), annotation, idx2goTerm,
                                         t2ICchilds.at(selectedGOs.at(t)));

                double tpD = tp / t2TPsum.at(selectedGOs.at(t));
                TPMatrix.at(t, pos) = tpD;
            }
        }
    }

    return TPMatrix;
}
*/

arma::Mat<double>
computeTPMatrix(Annotation &annotation,
                robin_hood::unordered_map<NetworKit::node, GOTerm> &idx2goTerm,
                std::vector<NetworKit::node> &selectedGOs) {

    arma::Mat<double> TPMatrix(selectedGOs.size(), selectedGOs.size(),
                               arma::fill::zeros);
    size_t t;
    robin_hood::unordered_map<NetworKit::node, double> t2TPsum;
    robin_hood::unordered_map<NetworKit::node, double> t2ICchilds;

    for (t = 0; t < selectedGOs.size(); t++) {
        size_t children_t =0;
        arma::Col<double> childIC(idx2goTerm.at(selectedGOs.at(t)).get_childrens().size());
        arma::Col<double> annotDiff(idx2goTerm.at(selectedGOs.at(t)).get_childrens().size());

        for (NetworKit::node c : idx2goTerm.at(selectedGOs.at(t)).get_childrens()) {
            childIC(children_t) =  idx2goTerm.at(c).get_IC();
            //icsum= icsum+idx2goTerm.at(c).get_IC();
            annotDiff(children_t) =  static_cast<double>(annotation.getGenes(c).size()) /
                    static_cast<double>(annotation.getGenes(selectedGOs.at(t)).size());
            children_t++;

        }

        double icsum = arma::sum(childIC);
        t2TPsum.insert(
        {selectedGOs.at(t), arma::sum(annotDiff + childIC/icsum)});
        t2ICchilds.insert({selectedGOs.at(t), icsum});
    }

    for (t = 0; t < selectedGOs.size(); t++) {

        for (NetworKit::node cn :
             idx2goTerm.at(selectedGOs.at(t)).get_childrens()) {
            if (find(selectedGOs.begin(), selectedGOs.end(), cn) !=
                    selectedGOs.end()) {
                size_t pos = static_cast<size_t>(
                            distance(selectedGOs.begin(),
                                     find(selectedGOs.begin(), selectedGOs.end(), cn)));
                double tp =
                        transitionalProb(cn, selectedGOs.at(t), annotation, idx2goTerm,
                                         t2ICchilds.at(selectedGOs.at(t)));

                double tpD = tp / t2TPsum.at(selectedGOs.at(t));
                TPMatrix.at(t, pos) = tpD;
            }
        }
    }
    return TPMatrix;
}



double TO(const robin_hood::unordered_set<NetworKit::node> &setGO1,
          const robin_hood::unordered_set<NetworKit::node> &setGO2,
          const int &max) {

    uint64_t ctr = 0; // Keep track of intersects
    for (const auto &item : setGO2) {
        if (setGO1.find(item) != setGO1.end()) {
            ctr++;
        }
    }

    return static_cast<double>(ctr) / static_cast<double>(max);
}

arma::Mat<double> computeSS(Annotation &annotation,
                            std::vector<std::string> &geneNodes, int &threads,
                            double &threshold) {
    std::cerr << "computeSS\n";
    arma::Mat<double> ssMatrix(geneNodes.size(), geneNodes.size(),
                               arma::fill::zeros);
    // BS: set diagonal elements to 0. Marginally faster than doing it in a loop
    ssMatrix.eye();

    // BS: Node IDs from the annotation might be whatever, so we remap them in
    // range [0, N] where N is the maximum number of annotations. We need two
    // objects to keep track of that in the next loop
    uint64_t new_index = 0;
    robin_hood::unordered_map<uint64_t, uint64_t> node_remap;

    // BS: Setup a sparse matrix [genes, GO_terms] that tracks which genes
    // are annotated to which GO terms
    arma::SpMat<uint64_t> gene_go_mat(geneNodes.size(), annotation.n_goterms());

    for (uint64_t i = 0; i < geneNodes.size(); i++) {
        if (annotation.findGene(geneNodes.at(i))) {
            // BS: Caching the result of .getGOTerms to use for both the matrix that
            // holds the annotations and the size. Always const & for speed and
            // correctness
            const auto &go_annotations = annotation.getGOTerms(geneNodes.at(i));
            for (const auto &term : go_annotations) {
                auto ptr = node_remap.find(term); // cache pointer
                uint64_t j = 0;
                if (ptr == node_remap.end()) {
                    // if it's a new GO term we assign the next column as its ID
                    // sequentially
                    node_remap[term] = new_index;
                    j = new_index;
                    new_index++;
                } else {
                    // otherwise we have seen it before and assign that column
                    j = ptr->second;
                }
                gene_go_mat(i, j) = 1;
            }
        }
    }

    // Drop unused columns
    std::cerr << "Sp origin: [" << gene_go_mat.n_rows << "," << gene_go_mat.n_cols
              << "]\n";
    gene_go_mat.resize(gene_go_mat.n_rows, new_index);
    std::cerr << "Sp resize: [" << gene_go_mat.n_rows << "," << gene_go_mat.n_cols
              << "]\n";

    std::cerr << "isec\n";

    // BS: Fast computation of all intersects with some matrix math
    arma::SpMat<uint64_t> gene_isec = gene_go_mat * gene_go_mat.t();
    std::cerr << "ISec Mat density: "
              << static_cast<double>(gene_isec.n_nonzero) /
                 static_cast<double>(gene_isec.n_rows * gene_isec.n_rows)
              << '\n';

    // BS: Keep diagonal as dense vector for faster random access. Diagonal has
    // the total number of annotation for a gene,
    auto isec_diag = arma::Col<uint64_t>(gene_isec.diag());

    // BS: Loop through nonzero elements and get the jaccard distance. Can be
    // optimized by looking at only the triangular matrix, but it seems using
    // arma iterators is still faster
    for (auto it = gene_isec.begin(); it != gene_isec.end(); ++it) {
        uint64_t i = it.row();
        uint64_t j = it.col();
        double v = static_cast<double>(*it);
        double m = isec_diag(i) > isec_diag(j) ? isec_diag(i) : isec_diag(j);
        ssMatrix(i, j) = ssMatrix(j, i) = v / m;
    }

    ssMatrix.elem(find(ssMatrix < threshold)).zeros();
    return ssMatrix;
}

arma::Mat<double> normilizeMatrix(arma::Mat<double> &matrix) {
    arma::rowvec d2 = arma::sum(matrix, 0);

    size_t s = d2.size();
    size_t i;
    arma::Mat<double> nMatrix = matrix;
    for (i = 0; i < s; i++) {
        double x = 1. / sqrt(d2.at(i));
        for (size_t j = i; j < s; j++) {
            double y = 1. / sqrt(d2.at(j));
            nMatrix(i, j) = nMatrix(j, i) = nMatrix(i, j) * (x * y);
        }
    }
    return nMatrix;
}

void exportVector(std::string filename, std::vector<NetworKit::node> &v) {

    std::stringstream ss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i != 0)
            ss << ",";
        ss << v[i];
    }
    std::ofstream myfile;
    myfile.open(filename);
    myfile << ss.str();
    myfile.close();
}

void exportVector(std::string filename, std::vector<std::string> &v) {

    std::stringstream ss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i != 0)
            ss << ",";
        ss << v[i];
    }
    std::ofstream myfile;
    myfile.open(filename);
    myfile << ss.str();
    myfile.close();
}

void preComputeMatrix(int threads, double threshold_ss, std::string obofile,
                      std::string goafile, std::string networkfile,
                      std::string outFolder) {
    std::vector<std::string> geneNodes;
    Annotation annotation;
    std::vector<NetworKit::node> selectedGOs;
    /*
   * Create Matrix
   */
    // ex netFile = "/home/aaron/git/Umea/GNI_predictors/ext-data/edgelist.txt"
    auto start = high_resolution_clock::now();
    arma::Mat<double> ppi = createPPIfromFile(networkfile, geneNodes);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "To read the network file took " << duration.count()
              << " seconds." << std::endl;
    // colvec c = sum(ppi,1);
    // cout<<"vector sum columns, position 0 = "<< c.at(0)<<endl;
    // READ GO file
    NetworKit::Graph go(
                0, true,
                true); // the first true is used to consider the graph as weigthed and the
    // second true to consider the graph directed
    robin_hood::unordered_map<NetworKit::node, GOTerm> idx2goterm;
    robin_hood::unordered_map<std::string, NetworKit::node> gotermsID2idx;
    robin_hood::unordered_map<std::string, Root> namespace2root;
    start = high_resolution_clock::now();
    // ex obofile = "/home/aaron/git/Umea/GNI_predictors/data/go.obo"
    ::reader(obofile, idx2goterm, go, gotermsID2idx, namespace2root, threads);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "To read the Gene Ontology took " << duration.count()
              << " seconds." << std::endl;
    // READ GOA file
    start = high_resolution_clock::now();
    // ex annotFile =
    // "/home/aaron/git/Umea/GNI_predictors/ext-data/Potra_gene_go.tsv"
    readAnnotation(goafile, annotation, gotermsID2idx, idx2goterm, selectedGOs);

    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "To read the Annotation file took " << duration.count()
              << " seconds." << std::endl;
    // Create association Matrix
    start = high_resolution_clock::now();
    arma::Mat<double> selAssoMatrix =
            computeGxTMatrix(geneNodes, selectedGOs, annotation);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "To create the gene2GOterm association matrix took "
              << duration.count() << " seconds." << std::endl;
    // Create the transitional probability matrix
    start = high_resolution_clock::now();
    arma::Mat<double> TPMatrix =
            computeTPMatrix(annotation, idx2goterm, selectedGOs);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "To create the transitional probability matrix took "
              << duration.count() << " seconds." << std::endl;
    // Create the semantic ppi
    start = high_resolution_clock::now();

    arma::Mat<double> ssMatrix =
            computeSS(annotation, geneNodes, threads, threshold_ss);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "To create the semantic ppi took " << duration.count()
              << " seconds." << std::endl;
    //std::cout << "TESTF" << std::endl;
    /*
   * NORMALIZE
   */
    //PPI
    start = high_resolution_clock::now();
    arma::Mat<double> nPpi = normilizeMatrix(ppi);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout <<"To normalize the PPI matrix took "<< duration.count()<<" seconds."
             << std::endl;
    // Semantic PPI
    start = high_resolution_clock::now();
    arma::Mat<double> nSSPpi = normilizeMatrix(ssMatrix);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout <<"To normalize the semantic PPI matrix took "<< duration.count()<<" seconds." << std::endl;
    start = high_resolution_clock::now();
    arma::Mat<double> hybrid = nPpi + nSSPpi;
    std::cout <<"Exporting..."
             << std::endl;
    /*
   * Save Matrix
   */
    hybrid.save( outFolder+"/hybrid.bin");
    nPpi.save( outFolder+"/nPpi.bin");
    nSSPpi.save( outFolder+"/nSSPpi.bin");
    TPMatrix.save(outFolder+"/TPMatrix.bin");
    selAssoMatrix.save(outFolder+"/annotationMatrix.bin");
 /*
    * Save vectors
    */
    exportVector(outFolder+"/genes.txt",geneNodes);
    exportVector(outFolder+"/GOterms.txt",selectedGOs);
    std::cout <<"Done."
}

int matrixPreparation(int ac, char *av[]) {
    int threads;
    double threshold_ss = 0.6;
    std::string outFolder;
    std::string oboFile;
    std::string goaFile;
    std::string networkFile;

    po::options_description desc("MatrixPreparation options");
    desc.add_options()("help,h", "produce help message")(
                "num_threads,T", po::value<int>(&threads)->default_value(1),
                "set number of threads to parallize")(
                "ss_threshold,S",
                po::value<double>(&threshold_ss)->default_value(static_cast<double>(0.6)),
                "Set a semantic similarity value threshold to filter links between "
                "genes/proteins with low similarity")(
                "goFile,G", po::value<std::string>(&oboFile),
                "Set the path of gene ontology file [OBO format]")(
                "annotationFile,A", po::value<std::string>(&goaFile),
                "Set the path of annotation File [tsv format]")(
                "networkInput,N", po::value<std::string>(&networkFile), "network file")(
                "outFolder,O",
                po::value<std::string>(&outFolder)->default_value("./results"),
                "output folder path to provide the results");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    po::store(po::command_line_parser(ac, av).options(desc).run(), vm);
    //  cout<<ac<<endl;
    if (vm.count("help") || ac == 2) {
        std::cerr << desc << '\n';
        return 1;
    }
    std::vector<std::string> err;
    bool flag = false;
    if (!vm.count("goFile")) {
        err.push_back("[ERROR] Please set the gene ontology file [OBO format] \n");
        flag = true;
    }

    if (!vm.count("annotationFile")) {
        err.push_back("[ERROR] Please set the annotation File [tsv format] \n");
        flag = true;
    }
    if (!vm.count("networkInput")) {
        err.push_back("[ERROR] Please set the network File\n");
        flag = true;
    }

    if (flag) {
        for (std::string s : err) {
            std::cerr << s;
        }
        std::cerr << desc << std::endl;
        return 1;
    }

    if (!boost::filesystem::exists(oboFile)) {
        err.push_back("[ERROR] Please use a correct path including the obo file.");
        flag = true;
    }
    if (!boost::filesystem::exists(goaFile)) {
        err.push_back("[ERROR] Please use a correct path including the goa file.");
        flag = true;
    }
    if (!boost::filesystem::exists(networkFile)) {
        err.push_back(
                    "[ERROR] Please use a correct path including the network file.");
        flag = true;
    }

    if (flag) {
        for (std::string s : err) {
            std::cerr << s;
        }
        std::cerr << desc << std::endl;
        return 1;
    }

    if (!boost::filesystem::exists(outFolder)) {
        boost::filesystem::create_directory(outFolder);
    }

    preComputeMatrix(threads, threshold_ss, oboFile, goaFile, networkFile,
                     outFolder);

    return 0;
}
