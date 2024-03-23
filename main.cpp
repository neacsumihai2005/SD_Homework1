/*
    baza 10
    (numere naturale)
*/
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

ifstream fin ("test.in");
ofstream fout ("test.out");

class Heap{
    ///minHeap
    /*
        Observatie:
        la Heap, in intervalul de pozitii
        [dim/2 + 1, dim]
        am numai frunze (indexare de la 1)
    */

private:
    int dim;
    int *v;

public:
    ~Heap(){
        delete v;
    }
    Heap(int N, int w[]){
        dim = N;
        v = new int[dim + 1];
        for(int i = 1; i <= N; i++){
            v[i] = w[i];
        }
        for(int i = dim / 2; i >= 1; i--){
            sift(i);
        }
    }

    int getDim(){
        return dim;
    }

    inline int left_son(int node){
        return node * 2;
    }
    inline int right_son(int node){
        return node * 2 + 1;
    }
    inline int father(int node){
        return node / 2;
    }
    inline int mini(){
        return v[1];
    }

    void sift(int node){
        /*
            ambii subarbori ai nodului au structura de heap corecta.
            eu "cern" nodul pana la locul potrivit, interschimband
            mereu cu cel mai mic fiu al sau
            pana ajunge sa fie mai mic decat ambii fii
            sau pana nu mai are fii xD
        */

        if(left_son(node) <= dim && right_son(node) <= dim){
            ///are ambii fii

            if(v[node] <= v[left_son(node)] && v[node] <= v[right_son(node)]){
                ///e mai mic decat ambii; ma opresc
                return;
            }
            else {
                if(v[left_son(node)] <= v[right_son(node)]){
                    swap(v[node], v[left_son(node)]);
                    sift(left_son(node));
                    return;
                }
                else {
                    swap(v[node], v[right_son(node)]);
                    sift(right_son(node));
                    return;
                }
            }
        }
        else if(left_son(node) <= dim){
            ///are doar un fiu, pe stanga

            if(v[node] <= v[left_son(node)]){
                ///e mai mic decat fiul; ma opresc
                return;
            }
            else {
                swap(v[node], v[left_son(node)]);
                sift(left_son(node));
                return;
            }
        }
        else {
            ///nu mai are fii
            return;
        }

    }

    void pop(){
        if(!(dim >= 1)){
            return;
        }

        swap(v[1], v[dim]);
        dim--;
        sift(1);
    }

    friend istream& operator >> (istream& in, Heap& X){
        in >> X.dim;
        X.v = new int[X.dim + 1];

        for(int i = 1; i <= X.dim; i++){
            in >> X.v[i];
        }
        for(int i = X.dim / 2; i >= 1; i--){
            X.sift(i);
        }



        return in;
    }

    friend ostream& operator << (ostream& out, Heap const& X){
        out << "Heap-ul are " << X.dim << " elemente: " << "\n";
        for(int i = 1; i <= X.dim; i++){
            out << X.v[i] << " ";
        }
        out << "\n";
    }
};

class Vector{
private:
    int dim;
    int *v;

public:
    friend class Heap;

    Vector(){}
    Vector(int x){
        dim = x;
        v = new int[dim + 1];
        for(int i = 0; i <= dim; i++){
            v[i] = 0;
        }
    }

    friend istream& operator >> (istream& in, Vector &X){
        in >> X.dim;
        X.v = new int[X.dim + 1];


        for(int i = 1; i <= X.dim; i++){
            in >> X.v[i];
        }

        return in;
    }

    friend ostream& operator << (ostream& out, Vector const &X){
        out << "Vectorul are " << X.dim << " elemente:" << "\n";
        for(int i = 1; i <= X.dim; i++){
            out << X.v[i] << " ";
        }
        out << "\n";
    }

    int getMax(){
        int mx = v[dim];
        for(int i = 1; i <= dim; i++){
            if(v[i] > mx){
                mx = v[i];
            }
        }
        return mx;
    }

    void countSort(int exp){
        int output[dim + 1];
        int ct[10] = {0};

        for(int i = 1; i <= dim; i++){
            int cifra = (v[i] / exp) % 10;
            ct[cifra]++;
        }

        for(int i = 1; i <= 9; i++){
            ct[i] = ct[i] + ct[i - 1];
        }

        for(int i = dim; i >= 1; i--){
            int cifra = (v[i] / exp) % 10;
            output[ ct[cifra] ]= v[i];
            ct[cifra]--;
        }

        for(int i = 1; i <= dim; i++){
            v[i] = output[i];
        }

    }

    void radixSort(){
        int mx = getMax();

        for(int exp = 1; exp <= mx; exp = exp * 10){
            countSort(exp);
            ///cout << *(this);

        }
    }


    ///functie care imi interclaseaza subvectorii:
    ///v1: subvectorul de la pozitia st la mij al vectorului tata
    ///v2: subvectorul de la pozitia mij + 1 la dr al vectorului tata
    ///in acelasi vector tata, pe intervalul [st, dr]
    void combina(int st, int mij, int dr){
        ///voi crea doi vectori temporari (indexati de la 1)
        int dimSt = mij - st + 1;
        int dimDr = dr - (mij + 1) + 1;
        auto *stArr = new int[dimSt + 1];
        auto *drArr = new int[dimDr + 1];

        for(int i = 1; i <= dimSt; i++){
            stArr[i] = v[st + i - 1];
        }
        for(int i = 1; i <= dimDr; i++){
            drArr[i] = v[(mij+1) + i - 1];
        }

        int indexSt = 1;
        int indexDr = 1;
        int indexV = st;
        while(indexSt <= dimSt && indexDr <= dimDr){
            if(stArr[indexSt] <= drArr[indexDr]){
                v[indexV] = stArr[indexSt];
                indexSt++;
            }
            else {
                v[indexV] = drArr[indexDr];
                indexDr++;
            }
            indexV++;
        }

        for(int i = indexSt; i <= dimSt; i++){
            v[indexV] = stArr[i];
            indexV++;
        }
        for(int i = indexDr; i <= dimDr; i++){
            v[indexV] = drArr[i];
            indexV++;
        }

        delete[] stArr;
        delete[] drArr;
    }

    void mergeSortUtil(int st, int dr){
        if(!(st < dr)){
            return;
        }

        int mij = st + (dr - st) / 2;
        mergeSortUtil(st, mij);
        mergeSortUtil(mij + 1, dr);
        combina(st, mij, dr);
    }
    void mergeSort(){
        mergeSortUtil(1, dim);
    }


    void shellSort(){
        for(int gap = dim / 2; gap >= 1; gap = gap / 2){
            for(int i = gap + 1; i <= dim; i++){
                int temp = v[i];

                int j;
                for(j = i; j - gap >= 1 && v[j - gap] > temp; j = j - gap){
                    v[j] = v[j - gap];
                }
                v[j] = temp;
            }
        }
    }

    void heapSortUtil(int st, int dr){
        int N = dr - st + 1;
        Heap X(N, v + st - 1);
        for(int i = st; i <= dr; i++){
            v[i] = X.mini();
            X.pop();
        }
    }
    void heapSort(){
        heapSortUtil(1, dim);
    }

    int partitie(int st, int dr){
        int pivot = v[dr];
        int lastBun = st - 1;
        for(int i = st; i <= dr; i++){
            if(v[i] < pivot){
                lastBun++;
                swap(v[i], v[lastBun]);
            }
        }
        swap(v[lastBun + 1], v[dr]);
        return lastBun + 1;
    }

    void quickSortUtil(int st, int dr){
        if(!(st < dr)){
            return;
        }

        int pi = partitie(st, dr);

        quickSortUtil(st, pi - 1);
        quickSortUtil(pi + 1, dr);
    }

    void quickSort(){
        quickSortUtil(1, dim);
    }

    void insertionSortUtil(int st, int dr){
        for(int i = st; i <= dr; i++){
            int temp = v[i];
            int k = i - 1;
            while(temp < v[k] && k >= st){
                v[k + 1] = v[k];
                k--;
            }
            v[k + 1] = temp;
        }
    }
    void insertionSort(){
        insertionSortUtil(1, dim);
    }

    void introSortUtil(int st, int dr, int depthLimit){
        if(!(st < dr)){
            return;
        }

        int N = dr - st + 1;

        if(N <= 16){
            insertionSortUtil(st, dr);
        }
        if(depthLimit == 0){
            heapSortUtil(st, dr);
        }
        else {
            int pi = partitie(st, dr);
            introSortUtil(st, pi - 1, depthLimit - 1);
            introSortUtil(pi + 1, dr, depthLimit - 1);
        }
    }

    void introSort(){
        int depthLimit = 2 * floor(log(dim));
        introSortUtil(1, dim, depthLimit);
    }
};



int main() {
    Vector X;
    fin >> X;

    ///X.mergeSort();
    ///X.radixSort();
    ///X.shellSort();
    ///X.heapSort();
    ///X.quickSort();
    ///X.insertionSort();
    X.introSort();

    fout << X;
    return 0;
}