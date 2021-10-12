#include "efficiency.h"
#include "makeHists.h"
#include <iostream>
#include <fstream>
#include "TH1D.h"

using namespace std;
int main(int argc, char **argv)
{
    //usage output
    if (argc != 4)
        std::cout << "usage : ./efficiency.exe filelist type filenumber" << std::endl;

    //read in filelist
    ifstream infile;
    infile.open(argv[1], ios::in);

    int count = 0;

    cout << "**********************************" << endl;
    cout << "*********efficiency Begin*********" << endl;
    cout << "**********************************" << endl;
    cout << endl;

    efficiency Run(argv[1], argv[2], argv[3]);

    TString inputRoot;
    while (infile >> inputRoot)
    {
        count++;
        Run.Initial(inputRoot, count, argv[2]);
        Run.Loop(argv[2]);
        Run.End(count);
    }

    infile.close();
    Run.Save(count);
    cout << endl;
    cout << "**********************************" << endl;
    cout << "************efficiency End************" << endl;
    cout << "**********************************" << endl;

    return 1;
}
