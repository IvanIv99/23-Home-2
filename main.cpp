#include <iostream>
#include "Matrix.h"

using namespace std;

int main() {
    cout << boolalpha;
    Matrix<float> Matrix1(2, 2);
    Matrix1.initialize();
    Matrix<float> Matrix2(1, 1);
    Matrix2.initialize();

    Matrix1 + Matrix2;
    Matrix1.print();

    cout << Matrix1.getDeterminant() << endl;
    !Matrix1;
    Matrix1.print();
    return 0;
}