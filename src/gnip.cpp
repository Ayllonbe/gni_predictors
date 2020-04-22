// External
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <matrixPreparation.h>
#include <randomWalk.h>
using namespace std;

const string usage_msg =
  "Gene function prediction.\n\n"
  "[Prepare matrix]\n"
  "  matrixPreparation         \t Create matrix to compute the prediction.\n"
  "\n"
  "[Run random Walk]\n"
  "  randomWalk                \t Run the random walk to get gene predictions.\n"
  "Beta version \n";


int main(int argc, char* argv[])
{
    try
      {
        if (argc < 2)
        {
          cerr << usage_msg << endl;
          throw std::invalid_argument("Too few arguments.");
        }
      }
      catch (std::invalid_argument& e)
      {
        cout << e.what() << '\n';
        return EPERM;
      }

    std::string task = argv[1];
    int ret;
    try
      {
        if (task == "matrixPreparation")
        {
          ret = matrixPreparation(argc, argv);
        }
        else if (task == "randomWalk")
        {
          ret =  newGOA(argc, argv);
        }
        else
        {
         cout << usage_msg << '\n';
          throw std::invalid_argument("Unrecognized task.");
        }
      }
      catch (std::exception& e)
      {
        cout << e.what() << '\n';
        return 1;
      }
      return ret;


}
