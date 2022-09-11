#ifndef SIMPLEX_H_INCLUDED
#define SIMPLEX_H_INCLUDED
#include<vector>
#include<fstream>
using namespace std;
class Rac
{
    int sz;
    int n;
    void simplify();
    public:
        Rac(int sz, int n);
        int lnko(int x, int y);
        int getsz()const;
        int getn()const;
        void setsz(int a);
        void setn(int a);
        void print()const;
        const Rac operator+(const Rac& x)const;
        const Rac& operator+=(const Rac& x);
        const bool operator==(const Rac& x) const;
        const Rac operator-(const Rac& x)const;
        const bool operator!=(const Rac& x)const;
        operator double()const;
        Rac reciprok();
        const Rac operator*(const Rac& x)const;
        const Rac operator*=(const Rac& x);
        const Rac operator/(Rac& x)const;
};
class Simplex
{
    vector<vector<Rac>> A;
    vector<vector<Rac>> PivotT;
    vector<Rac> b;
    vector<Rac> Z;
    vector<int> egyenlo;
    vector<int> nagyobb;
    vector <int> mmek;
    bool minn;
    double ertek;
    void korlell();
    void PivotTk();
    void printpivot();
    const unsigned int Pivotelem(const int x)const;
    void pivotalas(const unsigned int x,const unsigned int y);
    public:
        Simplex();
        ~Simplex();
        void filebol(const std::string& hely);
        void insert_cond(const std::string& szov);
        void insert_Z(const std::string& szov);
        void Aprint();
        void bprint();
        void Zprint();
        void printmeg();
        void solve();
        const double geter();
};
#endif // SIMPLEX_H_INCLUDED
