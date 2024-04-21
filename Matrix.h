#ifndef DZ2_MATRIX_H
#define DZ2_MATRIX_H

#include <functional>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

template<typename T>
class Matrix {
private:
    long int rows;
    long int cols;
    T** arr;

public:
    Matrix() : rows(0), cols(0), arr(nullptr) {}
    Matrix(long int m, long int n) {
        rows = m;
        cols = n;
        arr = new T*[rows];
        for (unsigned int i = 0; i < rows; ++i) {
            arr[i] = new T[cols];
            for (unsigned int j = 0; j < cols; ++j) {
                arr[i][j] = 0;
            }
        }
    }
    explicit Matrix(const char *filename) : rows(0), cols(0), arr(nullptr) {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Error");
        }
        std::string line;
        int i = 0;
        int rows_tmp = 0;
        int cols_tmp = 0;
        while (std::getline(file, line)) {
            if (i == 0) {
                std::stringstream ss(line);
                ss >> rows_tmp;
                ss.ignore(1, ',');
                ss >> cols_tmp;
                if (cols_tmp > 0 && rows_tmp > 0) {
                    rows = rows_tmp;
                    cols = cols_tmp;
                    arr = new T*[rows];
                    for (unsigned int l = 0; l < rows; ++l) {
                        arr[l] = new T[cols];
                    }
                } else {
                    throw std::out_of_range("Error");
                }
            } else {
                std::stringstream ss(line);
                for (unsigned int j = 0; j < cols_tmp; j++) {
                    T value;
                    ss >> value;
                    ss.ignore(1, ',');
                    arr[i-1][j] = value;
                }
            }
            i += 1;
        }
        file.close();
    }
    Matrix(const Matrix &second) {
        rows = second.rows;
        cols = second.cols;
        arr = new T*[rows];
        for (unsigned int i = 0; i < rows; ++i) {
            arr[i] = new T[cols];
            if (arr[i]) {
                for (unsigned int j = 0; j < cols; ++j) {
                    arr[i][j] = second.arr[i][j];
                }
            } else {
                throw std::runtime_error("Error");
            }
        }
    }
    Matrix(Matrix &&second) noexcept {
        rows = second.rows;
        cols = second.cols;
        arr = second.arr;
        second.arr = nullptr;
    }
    ~Matrix() {
        for (unsigned int i = 0; i < rows; ++i) {
            delete[] arr[i];
        }
        delete[] arr;
    }

    void initialize() {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                std::cin >> arr[i][j];
            }
        }
    }

    static Matrix<T> one(long int n) {
        Matrix<T> result(n,n);
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                if (i == j) {
                    result.arr[i][j] = 1;
                } else {
                    result.arr[i][j] = 0;
                }
            }
        }
        return result;
    }
    static Matrix<T> zero(long int m, long int n) {
        Matrix<T> result(m,n);
        for (unsigned int i = 0; i < m; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                result.arr[i][j] = 0;
            }
        }
        return result;
    }


    T getElement(long int i, long int j) const {
        if (i <= rows && j <= cols && i >= 0 && j >= 0) {
            return arr[i][j];
        } else {
            throw std::out_of_range("Error");
        }
    }

    T getDeterminant() {
        T sum = 0;
        if (rows == cols) {
            if (rows == 1) {
                return arr[0][0];
            } else if (rows == 2) {
                return arr[0][0] * arr[1][1] - arr[1][0] * arr[0][1];
            } else {
                for (unsigned int j = 0; j < cols; j++) {
                    sum += (j % 2 == 0 ? 1 : -1) * arr[0][j] * makeMinor(0, j).getDeterminant();
                }
            }
        } else {
            throw std::runtime_error("Error");
        }
        return sum;
    }

    void print() const {
        if (rows > 0 && cols > 0) {
            for (unsigned int i = 0; i < rows; ++i) {
                std::cout << "[";
                for (unsigned int j = 0; j < cols; ++j) {
                    std::cout << getElement(i, j) << " ";
                }
                std::cout << "]\n";
            }
        }
    }

    void readFile() {
        Matrix<T> tmp("text.txt");
        *this = tmp;
    }
    void writeFile() {
        std::ofstream file("text.txt");
        if (file.is_open()) {
            file << rows << "," << cols << "," << std::endl;
            for (unsigned int i = 0; i < rows; i++) {
                for (unsigned int j = 0; j < cols; j++) {
                    file << arr[i][j] << ",";
                }
                file << std::endl;
            }
            file.close();
        } else {
            throw std::runtime_error("Error");
        }
    }

    void operator=(const Matrix &second) {
        if (this != &second && this->arr != nullptr) {
            for (unsigned int i = 0; i < rows; ++i)
                delete[] arr[i];
            delete[] arr;
            rows = second.rows;
            cols = second.cols;
            arr = new T *[rows];
            for (unsigned int i = 0; i < rows; ++i) {
                arr[i] = new T[cols];
                if (arr[i]) {
                    for (unsigned int j = 0; j < cols; ++j) {
                        arr[i][j] = second.arr[i][j];
                    }
                } else {
                    throw std::runtime_error("Error");
                }
            }
        }
    }
    void operator+(const Matrix &second) const {
        if (rows == second.rows and cols == second.cols) {
            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = 0; j < cols; ++j) {
                    arr[i][j] = arr[i][j] + second.arr[i][j];
                }
            }
        } else {
            std::cerr << "Error\n";
        }
    }
    void operator-(const Matrix &second) const {
        if (rows == second.rows and cols == second.cols) {
            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = 0; j < cols; ++j) {
                    arr[i][j] = arr[i][j] - second.arr[i][j];
                }
            }
        } else {
            std::cerr << "Error\n";
        }
    }
    void operator*(const Matrix &second) const {
        double sum;
        if (cols == second.rows) {
            for (unsigned int i = 0; i < rows; i++) {
                for (unsigned int j = 0; j < second.cols; ++j) {
                    sum = 0;
                    for (unsigned int k = 0; k < cols; ++k) {
                        sum += arr[i][k] * second.arr[k][j];
                    }
                    arr[i][j] = sum;
                }
            }
        } else {
            std::cerr << "Error\n";
        }
    }
    void operator*(double factor) const {
        if (factor != 1.0) {
            for (unsigned int i = 0; i < rows; i++) {
                for (unsigned int j = 0; j < cols; j++) {
                    arr[i][j] *= factor;
                }
            }
        }
    }
    void operator*=(const Matrix &second) const { *this * second; }
    void operator*=(double factor) const { if (factor != 1.0) *this * factor; }
    void operator+=(const Matrix &second) const { *this + second; }
    bool operator==(const Matrix &second) {
        bool sign = true;
        if (rows == second.rows && cols == second.cols) {
            for (unsigned int i = 0; i < rows; i++) {
                for (unsigned int j = 0; j < cols; j++) {
                    sign *= (arr[i][j] == second.arr[i][j]);
                }
            }
            return sign;
        }
        return false;
    }
    bool operator!=(const Matrix &second) { return !(*this == second); }
    bool operator==(double scalyar) {
        bool sign = true;
        double num = arr[0][0];
        if (rows == cols) {
            for (unsigned int i = 0; i < rows; i++) {
                for (unsigned int j = 0; j < cols; j++) {
                    if (i == j)
                        sign *= (arr[i][j] == num && arr[i][j] == scalyar);
                    else {
                        sign *= (arr[i][j] == 0);
                    }
                }
            }
        }
        return sign;
    }
    bool operator!=(double scalyar) { return !(*this == scalyar); }

    void operator!() {
        Matrix Inverse(rows, cols);
        double det = this->getDeterminant();
        if (det != 0.00) {
            Inverse = this->makeAdjoint();
            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = 0; j < cols; ++j) {
                    this->arr[i][j] = Inverse.arr[i][j] / det;
                }
            }
        }
        else {
            std::cerr << "Error\n";
        }
    }


private:
    void SwapElements(long int i, long int j) {
        double temp = arr[i][j];
        arr[i][j] = arr[j][i];
        arr[j][i] = temp;
    }
    void swapLines(long int i1, long int i2) {
        double *temp = arr[i1];
        arr[i1] = arr[i2];
        arr[i2] = temp;
    }

    void dif2(const double *a) {
        unsigned int k = 0;
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                arr[i][j] = a[k];
                k += 1;
            }
        }
    }

    void transp() {
        if (rows == cols) {
            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = i; j < cols; ++j) {
                    if (j != i)
                        this->SwapElements(i, j);
                }
            }
        } else {
            std::cerr << "Error\n";
        }
    }

    Matrix makeAdjoint() {
        Matrix Adjoint(rows, cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                Adjoint.arr[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * makeMinor(i, j).getDeterminant();
            }
        }
        Adjoint.transp();
        return Adjoint;
    }
    Matrix makeMinor(long int i, long int j) {
        int k = 0;
        Matrix Minor(rows - 1, cols - 1);
        double arr1[Minor.rows * Minor.cols];
        for (unsigned int brr_i = 0; brr_i <= Minor.rows; ++brr_i) {
            for (unsigned int brr_j = 0; brr_j <= Minor.cols; ++brr_j) {
                if (brr_i != i && brr_j != j) {
                    arr1[k] = arr[brr_i][brr_j];
                    k += 1;
                }
            }
        }
        Minor.dif2(arr1);
        return Minor;
    }
};


#endif //DZ2_MATRIX_H
