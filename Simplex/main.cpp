#include <iostream>
#include "simplex.h"

int main()
{
    Simplex proba;
    proba.insert_cond("1 0.3 0.5 -1 3 0 >= 100");
    proba.insert_cond("5 3 2 0 0 2 = 650");
    proba.insert_cond("1 0.5 1 3 2 1 <= 150");
    proba.insert_Z("max 100 160 250 0 15 30");
    proba.solve();
    std::cout<<proba.geter();
    std::cout<<std::endl;
    Simplex uj;
    uj.filebol("./text.txt");
    std::cout<<std::endl;
    uj.Aprint();
    uj.Zprint();
    uj.bprint();
    uj.solve();
    std::cout<<uj.geter();
    return 0;
}
