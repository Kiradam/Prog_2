#include  "simplex.h"
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <fstream>

    void Rac::simplify(){
	    if(n<0){sz=-sz;n=-n;}
        int t=n;
        while(sz%t||n%t)--t;
        sz/=t;n/=t;
    }
        Rac::Rac(int sz=0, int n=1){
            this->sz=sz; this->n=n;
            if(n==0) throw std::domain_error("0 a nevezoben");
            simplify();
            if((*this).n<0){this->sz*=-1;this->n*=-1;}
        }
        int Rac::lnko(int x, int y){
            if (x==0){return 1;}
            else if(x%y==0){
                if (x>y){return y;}
                else{return x;}}
            else if(x>y){



                return lnko(y,x%y);
            }
            else if(y>x){
                return lnko(y,x);
            }
            return 1;
        }
        int Rac::getsz()const{return sz;}
        int Rac::getn ()const{return n;}
        void Rac::setsz(int a){sz=a;simplify();}
        void Rac::setn (int a){n=a;
        if(n==0) throw std::domain_error("0 a nevezoben");
            simplify();
        }
        void Rac::print()const{
            if (sz%n==0){
            cout<<sz/n;}
            else{
            cout<<sz<<"/"<<n;}
        }
        const Rac Rac::operator+(const Rac& x)const {
            return Rac(sz*x.n+n*x.sz,n*x.n);
        }
        const Rac& Rac::operator+=(const Rac& x) {
            *this=*this+x;
            return *this;
        }
        const bool Rac::operator==(const Rac& x) const{
            if (sz==x.sz&&n==x.n){return true;}
            else{return false;}
        }
        const Rac Rac::operator-(const Rac& x)const{
            return Rac(sz*x.n-n*x.sz,n*x.n);
        }
        const bool Rac::operator!=(const Rac& x)const{
            if ((*this)==x){return false;}
            else {return true;}
        }
        Rac::operator double()const{return (double)sz/(double)n;}
        Rac Rac::reciprok(){
            if (this->sz==0){throw std::domain_error("0 a nevezoben");}
            return Rac(n,sz);}
        const Rac Rac::operator*(const Rac& x)const{
            return Rac(sz*x.sz,n*x.n);
        }
        const Rac Rac::operator*=(const Rac& x){
            *this=*this*x;
            return *this;
        }
        const Rac Rac::operator/(Rac& x)const{
            return Rac(sz*x.n,n*x.sz);
        }

    void Simplex::printpivot(){
            unsigned int i=0;
            unsigned int j=0;
            while(i<PivotT.size()){
                j=0;
                while(j<PivotT[0].size()){
                    PivotT[i][j].print();
                    cout<<" ";
                    j++;}
                cout<<'\n';
                i++;
            }
        }
    void Simplex::korlell(){
        unsigned int i=0;
        unsigned int j=0;
        bool nemko=true;
        bool nincsmo=true;
        while(i<mmek.size()){
            nemko=true;
            nincsmo=true;
            j=0;
            if(PivotT[PivotT.size()-1][mmek[i]]>0){
                while(j<A.size()){
                    if(PivotT[j][mmek[i]]>0){nincsmo=false;nemko=false;}
                j++;
            }
            j=0;
            if(nincsmo){
                nemko=false;
                while(j<A.size()){
                    if(PivotT[j][mmek[i]].getsz()==0){nemko=true;}
                j++;
                }
                if (nemko){printpivot();throw std::domain_error("Nem korlatos celfuggveny");}
                else if (nincsmo){printpivot();throw std::domain_error("Ellentmondasos rendszer.");}
            }
            }
            i++;
        }}
    void Simplex::PivotTk(){
            if(Z.size()!=A[0].size()){throw std::domain_error("Celfuggveny egyutthatovektora, es dontesi vektor dimenzioja nem azonos.");}
            unsigned int i=0;
            unsigned int j=0;
            bool megengedett=true;
            while (i<A.size()){
                PivotT.push_back(A[i]);
                i++;
            }
            if(nagyobb.size()!=0){
                i=1;
                PivotT.push_back(A[nagyobb[0]-1]);
                while(i<nagyobb.size()){
                    j=0;
                    while(j<A[0].size()){
                        PivotT[PivotT.size()-1][j]+=PivotT[nagyobb[i]-1][j];
                        j++;
                    }
                    i++;
                }
            }
            if(egyenlo.size()!=0){
                if(nagyobb.size()!=0){i=0;}
                else{PivotT.push_back(A[egyenlo[0]-1]);
                i=1;}
                while(i<egyenlo.size()){
                    j=0;
                    while(j<A[0].size()){
                        PivotT[PivotT.size()-1][j]+=PivotT[egyenlo[i]-1][j];
                        j++;
                    }
                    i++;
                }
            }
            PivotT.push_back(Z);
            i=0;
            while(i<A.size()){//+ketf
                j=0;
                while(j<A.size()){
                    if (j==i){
                        PivotT[i].push_back(Rac(1,1));
                    }
                    else{PivotT[i].push_back(Rac(0,1));}
                    j++;
                }
                i++;
            }
            j=0;
            while(j<=(A.size()+nagyobb.size())){
                PivotT[PivotT.size()-1].push_back(Rac(0,1));
                if(nagyobb.size()!=0||egyenlo.size()!=0){
                    PivotT[PivotT.size()-2].push_back(Rac(0,1));}
                j++;
            }
            i=0;
            j=0;
            while(i<nagyobb.size()){
                j=0;
                while (j<A.size()){
                    if (nagyobb[i]-1u==j){
                        PivotT[j].push_back(Rac(-1,1));
                    }
                    else{PivotT[j].push_back(Rac(0,1));}
                    j++;
                }
                i++;
            }
            i=0;
            while(i<b.size()){
                PivotT[i].push_back(b[i]);
                i++;
            }
            i=0;
            if (nagyobb.size()!=0||egyenlo.size()!=0){
            while(i<nagyobb.size()){
                PivotT[PivotT.size()-2][PivotT[0].size()-1]+=b[nagyobb[i]-1];
                i++;
            }
            i=0;
            while(i<egyenlo.size()){
                PivotT[PivotT.size()-2][PivotT[0].size()-1]+=b[egyenlo[i]-1];
                i++;
            }}
            i=0;
            while(i<nagyobb.size()){
                PivotT[PivotT.size()-2][A[0].size()+A.size()+i]+=Rac(-1,1);
                i++;
            }
            i=0;
            j=0;
            while(i<PivotT[0].size()-1){
                    megengedett=true;
                        j=0;
                        while(j<nagyobb.size()){
                            if(A[0].size()+nagyobb[j]-1==i){megengedett=false;}
                            j++;
                        }
                        j=0;
                        while(j<egyenlo.size()){
                            if(A[0].size()+egyenlo[j]-1==i){megengedett=false;}
                            j++;
                        }
                    if( megengedett) {mmek.push_back(i);}
                    i++;
                }
    }
    const unsigned int Simplex::Pivotelem(const int x)const{
            vector<int> megengedett;
            bool korl=false;
            bool nemures=false;
            unsigned i=0;
            int eredmeny;
            while(i<A.size()){
                if(PivotT[i][x]>=0){
                    nemures=true;
                    if(PivotT[i][x].getsz()!=0){korl=true;megengedett.push_back(i);}
                }
                i++;
            }
            if(nemures==false){throw std::domain_error("Nem korlatos celfugggveny");}
            else if(korl==false){throw std::domain_error("Nem korlatos celfuggveny");}
            i=1;
            eredmeny=megengedett[0];
            while(i<megengedett.size()){
                if((PivotT[eredmeny][PivotT[0].size()-1]/PivotT[eredmeny][x])>PivotT[megengedett[i]][PivotT[0].size()-1]/PivotT[megengedett[i]][x]&&PivotT[megengedett[i]][PivotT[0].size()-1]/PivotT[megengedett[i]][x]>0){eredmeny=megengedett[i];}
                i++;
            }
            megengedett.clear();
            return eredmeny;
        }
        void Simplex::pivotalas(const unsigned int x,const unsigned int y){
            unsigned int i=0;
            unsigned int j=0;
            Rac konstans=PivotT[x][y];
            while(i<PivotT[0].size()){
                PivotT[x][i]=PivotT[x][i]/konstans;
                i++;
            }
            i=0;
            while(i<PivotT.size()){
                if(i!=x){
                    j=0;
                    konstans=PivotT[i][y]/PivotT[x][y];
                    while(j<PivotT[0].size()){
                        PivotT[i][j]=PivotT[i][j]-(PivotT[x][j]*konstans);
                        j++;
                    }
                }
                i++;
            }
    }
        Simplex::Simplex(){minn=false;}
        Simplex::~Simplex(){Z.clear();b.clear();nagyobb.clear();mmek.clear();egyenlo.clear();
            unsigned int i=0;
            while (i<A.size()){
                A[i].clear();
                i++;
            }
            while (i<PivotT.size()){
                PivotT[i].clear();
                i++;
            }
            A.clear();
            PivotT.clear();
        }
        void Simplex::filebol(const std::string& hely){
            ifstream in(hely);
            std::string str;
            while(std::getline(in,str)){
                if(str[0]!='m'){
                    Simplex::insert_cond(str);}
                else{
                    Simplex::insert_Z(str);
                }
            }
        }
        void Simplex::insert_cond(const std::string& szov =NULL){;
            int i=0;
            int j=0;
            double szam;
            int ht=0; //hany tizedesjegy
            int elsz=0;//elozo szokoz helye
            char* ideigl=NULL;
            vector<Rac> beszur;
            if (szov.size()!=0){
                while(true){
                    ht=0;
                    if (szov[i]==' '){
                        ideigl=new char[i-elsz+1];
                        while (j<i-elsz){
                            if(szov[elsz+1+j]=='.'){ht=i-j-elsz-2;}
                            ideigl[j]=szov[elsz+j];
                            j++;
                        }
                        ideigl[j]='\0';
                        if(ideigl[1]=='<'||ideigl[1]=='>'||ideigl[1]=='='){
                            if(ideigl[1]=='>'){nagyobb.push_back(A.size()+1);}
                            else if(ideigl[1]=='='){egyenlo.push_back(A.size()+1);}
                        }
                        else{sscanf(ideigl, "%lf", &szam);
                            beszur.push_back(Rac(int(szam*pow(10,ht)),int(pow(10,ht))));
                        }
                        delete [] ideigl;
                        j=0;
                        elsz=i;
                    }
                    i++;
                    if(szov[i]=='\0'){break;}
                }
                if(A.size()!=0 && A[0].size()!=beszur.size()){
                    throw std::domain_error("Nem megfelelo alaku a feltetel!");
                    }
                A.push_back(beszur);
                beszur.clear();
            }
            else{throw std::domain_error("Nem adtal meg beolvasando adatot!");}
            j=0;
            ideigl=new char[i-elsz+1];
            while (j<=i-elsz){
                    if(szov[elsz+1+j]=='.'){ht=i-j-elsz-2;}
                    ideigl[j]=szov[elsz+j];
                    j++;
            }
            sscanf(ideigl, "%lf", &szam);
            if (szam<0){
                throw std::domain_error("A jobb oldal negativ!!!");
            }
            b.push_back(Rac(int(szam*pow(10,ht)),int(pow(10,ht))));;
            delete [] ideigl;
        }
        void Simplex::insert_Z(const std::string& szov =NULL){;
            int i=4;
            int j=0;
            unsigned int k=0;
            double szam;
            int ht=0; //hany tizedesjegy
            int elsz=4;//elozo szokoz helye
            char* ideigl=NULL;
            vector<Rac> beszur;
            if (szov.size()!=0){
                while(true){
                    ht=0;
                    if (szov[i]==' '){
                        ideigl=new char[i-elsz+1];
                        while (j<i-elsz){
                            if(szov[elsz+1+j]=='.'){ht=i-j-elsz-2;}
                            ideigl[j]=szov[elsz+j];
                            j++;
                        }
                        ideigl[j]='\0';

                        sscanf(ideigl, "%lf", &szam);
                            beszur.push_back(Rac(int(szam*pow(10,ht)),int(pow(10,ht))));

                        delete [] ideigl;
                        j=0;
                        elsz=i;
                    }
                    i++;
                    if(szov[i]=='\0'){break;}
                }
                Z=beszur;
                beszur.clear();
            }
            j=0;
            ideigl=new char[i-elsz+1];
            while (j<=i-elsz){
                    if(szov[elsz+1+j]=='.'){ht=i-j-elsz-2;}
                    ideigl[j]=szov[elsz+j];
                    j++;
            }
            sscanf(ideigl, "%lf", &szam);
            Z.push_back(Rac(int(szam*pow(10,ht)),int(pow(10,ht))));
            if(szov[1]=='i'){
                minn=true;
                while(k<Z.size()){
                    Z[k]*=Rac(-1,1);
                    k++;
                }
            }
            delete [] ideigl;
        }
        void Simplex::Aprint(){
            unsigned int i=0;
            unsigned int j=0;
            while(i<A.size()){
                j=0;
                while(j<A[0].size()){
                A[i][j].print();
                cout<<" ";
                j++;}
                cout<<'\n';
                i++;
            }
        }
        void Simplex::bprint(){
            unsigned int i=0;
            while(i<b.size()){
                b[i].print();
                cout<<" ";
                i++;
            }
            cout<<endl;
        }
        void Simplex::Zprint(){
            unsigned int i=0;
            while(i<Z.size()){
                Z[i].print();
                cout<<" ";
                i++;
            }
            cout<<endl;
        }
        /*void printpivot(){
            unsigned int i=0;
            unsigned int j=0;
            while(i<PivotT.size()){
                j=0;
                while(j<PivotT[0].size()){
                    PivotT[i][j].print();
                    cout<<" ";
                    j++;}
                cout<<'\n';
                i++;
            }
        }*/
        void Simplex::printmeg(){
            unsigned i=0;
            while(i<mmek.size()){
            cout<<mmek[i]<<" ";
            i++;}
        }
        void Simplex::solve(){
            bool vanpoz;
            unsigned int i=0;
            unsigned int j=0;
            bool notsolved=true;
            if(PivotT.size()==0){PivotTk();}
            cout<<"indulasi tabla: "<<endl;
            printpivot();
            cout<<"Bazisba kerulhet: ";
            printmeg();
            cout<<endl;
            if(nagyobb.size()==0u&&egyenlo.size()==0u){
                while(notsolved){
                    i=0;
                    while(i<PivotT[0].size()){
                        if(PivotT[PivotT.size()-1][i].getsz()>0){
                            korlell();
                            pivotalas(Pivotelem(i),i);
                            break;
                        }
                        i++;
                    if (i==PivotT[0].size()-1){notsolved=false;}
                    }
                }
            }
            else{
                cout<<endl;
                while(notsolved){
                    i=0;
                    while(i<mmek.size()){
                        if(PivotT[PivotT.size()-2][mmek[i]].getsz()>0&&PivotT[PivotT.size()-2][PivotT[0].size()-1].getsz()!=0){
                            //printpivot();
                            //cout<<endl;
                            korlell();
                            pivotalas(Pivotelem(mmek[i]),mmek[i]);
                            //printpivot();
                            //cout<<endl;
                            break;
                        }
                        i++;
                    if ((PivotT[PivotT.size()-2][PivotT[0].size()-1].getsz()==0)){notsolved=false;}
                    else if(i==mmek.size()&&notsolved&&PivotT[PivotT.size()-2][PivotT[0].size()-1].getsz()!=0){printpivot();throw std::domain_error("Nincs kiindulo bazis");}
                    }
                }
                PivotT.erase(PivotT.end()-2);
                notsolved=true;
                while(notsolved){
                    i=0;
                    while(i<mmek.size()){
                        if(PivotT[PivotT.size()-1][mmek[i]].getsz()>0){
                            j=0;
                            vanpoz=false;
                            //korlell();
                            //printpivot();
                            //cout<<endl;
                            pivotalas(Pivotelem(mmek[i]),mmek[i]);
                            while(j<mmek.size()){
                                if(PivotT[PivotT.size()-1][mmek[j]].getsz()>0){vanpoz=true;}
                                j++;
                            }
                            break;
                        }
                        i++;
                    if (i==mmek.size()-1&&vanpoz==false){notsolved=false;}
                    }
                }
            }
            i=0;
            cout<<"Lellasi tabla: "<<endl;
            printpivot();
            ertek=PivotT[PivotT.size()-1][PivotT[0].size()-1]*Rac(-1,1);
        }
        const double Simplex::geter(){
            if (minn){return -ertek;}
            else {return ertek;}
        }

