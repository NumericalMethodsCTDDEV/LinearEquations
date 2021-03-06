#include <vector>
#include <cmath>
#include "linearSystemsSolver.h"

bool checkDiagPrev(std::vector < std::vector<double> > a,int amountOfVars,int amountOfEq){
    //проверка сходимости по методу - условие диагонального преобладания sum|aij|<|aii|

    for(int i=0;i<amountOfEq;i++) {
        double sum = 0;
        double diag=0;
        for (int j = 0; j < amountOfVars; j++)
            if (i==j)
                diag=std::fabs(a[i][j]);
            else
                sum+=std::fabs(a[i][j]);
        if (diag<=sum)
            return false;
    }
    return true;
}

void solveJS(bool isConv,std::vector < std::vector<double> > a,int amountOfVars,int amountOfEq,double* X,double* suppX, double EPS,std::vector<double> &ans){
    int counter =0;
    //изначальный вектор
    for (int j=0;j<amountOfVars;j++)
        X[j] = 0;

    double norm; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.
    ans.assign(amountOfVars,0);
    do {
        for (int i = 0; i < amountOfEq; i++) {
            suppX[i] = a[i][amountOfVars];//свободные члены
            for (int j = 0; j < amountOfVars; j++) {
                if (i != j)
                    suppX[i] -= a[i][j] * X[j];
            }
            suppX[i] /= a[i][i];
        }
        norm = std::fabs(ans[0] - suppX[0]);
        for (int j = 0; j < amountOfVars; j++) {
            if (std::fabs(ans[j] - suppX[j]) > norm)
                norm = std::fabs(ans[j] - suppX[j]);
            ans[j] = suppX[j];
        }
        std::swap(X,suppX);
        counter++;
    } while (norm > EPS&&(isConv||counter<linearSystemsSolver::ITERATIONS));
}

int jacobi(std::vector < std::vector<double> > a, std::vector<double> &ans)
{
    int amountOfEq = (int) a.size();
    int amountOfVars = (int) a[0].size() - 1;

    //проверка сходимости
    bool diagPrev=checkDiagPrev(a,amountOfVars,amountOfEq);

    //проверка корректности размеров матрицы
    /*if (amountOfEq!=amountOfVars)
        return 2;*/

    //вычисление апостериорной оценки
    std::vector < std::vector<double> > b(a);
    for(int i=0;i<amountOfEq;i++)
        for(int j=0;j<amountOfVars;j++)
            if (i==j)
                b[i][j]=0;
            else
                b[i][j]=(-1)*(a[i][j]/a[i][i]);

    double eps=linearSystemsSolver::norm(b);
    if (eps>1/2)
        eps=((1-eps)/eps)*linearSystemsSolver::EPS;
    else
        eps=linearSystemsSolver::EPS;

    //решение системы
    double* X = new double[amountOfVars];
    double* suppX = new double[amountOfVars];
    solveJS(diagPrev,a,amountOfVars,amountOfEq,X,suppX,eps,ans);
    delete[] suppX;
    delete[] X;

    if (diagPrev)
        return 1;
    else
        return 2;
}


/*int seidel(std::vector < std::vector<double> > a, std::vector<double> &ans)
{
    int amountOfEq = (int) a.size();
    int amountOfVars = (int) a[0].size() - 1;


    //проверка сходимости
    bool diagPrev=checkDiagPrev(a,amountOfVars,amountOfEq);

    //вычисление апостериорной оценки
    std::vector < std::vector<double> > b(a);//std::vector< std::vector<double>>(std::vector<double >(amountOfEq),amountOfVars);
    for(int i=0;i<amountOfEq;i++)
        for(int j=0;j<amountOfVars;j++)
            if (i==j)
                b[i][j]=0;
            else
                b[i][j]=(-1)*(a[i][j]/a[i][i]);

    double eps=linearSystemsSolver::norm(b);
    if (eps>1/2)
        eps=((1-eps)/eps)*linearSystemsSolver::EPS;
    else
        eps=linearSystemsSolver::EPS;

    //решение системы
    double* X = new double[amountOfVars];//тк сразу записывать в ответ
    solveJS(diagPrev,a,amountOfVars,amountOfEq,X,X,linearSystemsSolver::EPS,ans);
    delete[] X;
    if(diagPrev)
        return 1;
    else
        return 2;
}
*/
