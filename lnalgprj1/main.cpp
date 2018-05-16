#include<iostream>
#include<Windows.h>

using namespace std;

class Matrix {
private:
	double *matrix;
	double *identity;
	double *matvec;
	int row, col;

	void makeIdentity() {
		if (identity) delete identity;
		identity = new double[row*col]();
		for (int i = 0; i < row; i++) {
			identity[i*col + i] = 1;
		}
	}

	double& getIdmat(int row, int col) {
		if (identity != NULL && row >= this->row && col >= this->col)
			cout << "index out of range error" << endl;
		else
			return identity[row*this->col + col];
	}

	void globarmul(int row, double multiplier) {
		for (int i = 0; i < col; i++) {
			this->get(row, i) *= multiplier;
			if (identity)
				this->getIdmat(row, i) *= multiplier;
		}
		if (matvec)
			matvec[row] *= multiplier;
	}

	void changeline(double* mat, int destrow, int srcrow) {
		double tmp;
		for (int i = 0; i < col; i++) {
			tmp = mat[destrow*col+i];
			mat[destrow*col + i] = mat[srcrow*col + i];
			mat[srcrow*col + i] = tmp;
		}
	}

	void pivot(int row) {
		int idx = row;
		for (int i = row+1; i < this->row; i++) {
			if (this->get(i, row) != 0.0 && abs(this->get(i, row)) < abs(this->get(idx, row)))
				idx = i;
		}
		if (idx == row) return;
		changeline(this->matrix, idx, row);
		if (identity)
			changeline(this->identity, idx, row);
		if (matvec) {
			double tmp = matvec[idx];
			matvec[idx] = matvec[row];
			matvec[row] = tmp;
		}
	}

	double getmultiplier(int destrow, int srcrow) {
		return this->get(destrow, srcrow) / this->get(srcrow, srcrow);
	}

	void sub(int destrow, int srcrow, double multiplier) {
		for (int i = 0; i < col; i++) {
			this->get(destrow, i) -= this->get(srcrow, i)* multiplier;
			if (identity) 
				this->getIdmat(destrow, i) -= this->getIdmat(srcrow, i)* multiplier;
		}
		if (matvec)
			matvec[destrow] -= matvec[srcrow] * multiplier;
	}

	double determinant() {
		if (row != col) {
			cout << "row must be same column length" << endl;
			return 0;
		}
		if (row == 2) 
			return get(0, 0)*get(1, 1) - get(0, 1)*get(1, 0);
		double det = 0;
		for (int i = 0; i < col; i++) {
			Matrix *minor = minormatrix(0, i);
			if (i % 2 == 0)
				det += get(0, i)*minor->determinant();
			else
				det -= get(0, i)*minor->determinant();
			delete minor;
		}
		return det;
	}

	Matrix* minormatrix(int row, int col) {
		double *temp = new double[(this->row - 1)*(this->col - 1)];
		int counter = 0;
		for (int x = 0; x < this->row; x++) {
			for (int y = 0; y < this->col; y++) {
				if (x == row || y == col) continue;
				temp[counter] = this->get(x, y);
				counter++;
			}
		}
		Matrix *minor = new Matrix(this->row - 1, this->col - 1, temp);
		delete[] temp;
		return minor;
	}

public:
	Matrix(int row, int col) {
		matrix = new double[row*col]();
		this->row = row;
		this->col = col;
	}

	Matrix(int row, int col, const double mat[]) {
		matrix = new double[row*col];
		this->row = row;
		this->col = col;
		memcpy(matrix, mat, sizeof(double)*row*col);
	}

	Matrix(int row, int col, const double mat[], const double c_vec[]) : Matrix(row, col, mat) {
		matvec = new double[row];
		memcpy(matvec, c_vec, sizeof(double)*row);
	}

	void gauss() {
		for (int i = 0; i < row - 1; i++) {
			pivot(i);
			for (int x = i + 1; x < row; x++) {
				sub(x, i, getmultiplier(x, i));
			}
		}
	}

	void gaussjordan() {
		makeIdentity();
		gauss();
		//revert gauss
		for (int i = row - 1; i >= 0; i--) {
			for (int x = i - 1; x >= 0; x--) {
				sub(x, i, getmultiplier(x, i));
			}
			globarmul(i, 1.0 / get(i, i));
		}
	}

	Matrix* getIdentity() {
		if (!identity) return NULL;
		return new Matrix(row, col, identity);
	}

	int getRow() {
		return row;
	}

	int getCol() {
		return col;
	}

	double& get(int row, int col) { // it can modify actuall matrix!!
		if (row >= this->row && col >= this->col)
			cout << "index out of range error" << endl;
		else
			return matrix[row*this->col + col];
	}
	
	Matrix* operator*(Matrix& m) {
		if (this->col != m.row) {
			cout << "can not multiplication. check row, column size" << endl;
			return NULL;
		}
		Matrix *res = new Matrix(this->row, m.col);
		for (int x = 0; x < row; x++) {
			for (int y = 0; y < m.col; y++) {
				for(int i=0; i < col; i++){
					res->get(x, y) += this->get(x, i) * m.get(i, y);
				}
			}
		}
		return res;
	}
	
	friend ostream& operator<<(ostream& os, Matrix& mat) {
		os.precision(3);
		for (int x = 0; x < mat.row; x++) {
			for (int y = 0; y < mat.col; y++) {
				os << mat.get(x, y) << "\t";
			}
			if (mat.identity) {
				cout << "|\t";
				for (int y = 0; y < mat.col; y++) {
					os << mat.getIdmat(x, y) << "\t";
				}
			}
			if (mat.matvec) {
				cout << "|\t" << mat.matvec[x];
			}
			os << "\n";
		}
		return os;
	}

	void transpose() {
		double *replace = new double[row*col];
		for (int x = 0; x < row; x++) {
			for (int y = 0; y < col; y++) {
				replace[y*row+x] = get(x, y);
			}
		}
		delete[] matrix;
		matrix = replace;
		int temp = row;
		row = col;
		col = temp;
	}

	Matrix* getVector() {
		if (matvec == NULL) return NULL;
		return new Matrix(row, 1, matvec);
	}

	Matrix* getInverse() {
		double det = determinant();
		if (det == 0.0) return NULL;
		Matrix *res = new Matrix(row, col);
		for (int x = 0; x < row; x++) {
			for (int y = 0; y < col; y++) {
				Matrix *minor = minormatrix(x, y);
				if ((x + y) % 2 == 0)
					res->get(x, y) = minor->determinant() / det;
				else
					res->get(x, y) -= minor->determinant() / det;
				delete minor;
			}
		}
		res->transpose();
		return res;
	}

	double getDeterminant() {
		return determinant();
	}
	
	~Matrix() {
		if (matrix) delete[] matrix;
		if (identity) delete[] identity;
		if (matvec) delete[] matvec;
	}
};

int main(int argc, char *argv[]) {
	double A[][25] = { {
		1,2,3,1,2,
		2,1,3,1,3,
		3,5,2,4,2,
		2,3,1,3,1,
		3,4,5,2,1
	},
	{
		1,2,3,1,2,
		2,1,3,1,3,
		3,3,6,2,5,
		2,3,1,3,1,
		3,4,5,2,1
	},
	{
		1,2,3,1,2,
		2,1,3,1,3,
		3,3,6,2,5,
		2,3,1,3,1,
		5,4,5,2,1
	} };
	double Ax[][5] = {
		{ -4, -4, 1, 0, -7 },
		{ -4, -4, -8, 0, -7 },
		{ -4, -4, 6, 0, 7 }
	};

	for (int i = 0; i < 3; i++) {
		Matrix a(5, 5, A[i], Ax[i]);
		double det = a.getDeterminant();
		cout << "orignal matrix" << endl;
		if (det == 0) { 
			cout << a << endl;;
			cout << "determinant is 0. skip matrix." << endl;
			continue; 
		}
		Matrix *inv = a.getInverse();
		Matrix c(5, 1, Ax[i]);
		cout << a << endl;
		cout << "1. gause jordan (with inv)" << endl;
		a.gaussjordan();
		cout << a << endl;
		cout << "2. inverse matrix (with det)" << endl;
		cout << *inv << endl;
		Matrix *x = (*inv)*c;
		cout << "inv(A)*c" << endl;
		cout << *x << endl;
		cout << "3.Cramer" << endl;
		cout << "det(A) = " << det << endl;
		for (int x = 0; x < 5; x++) {
			Matrix a(5, 5, A[i], Ax[i]);
			for (int y = 0; y < 5; y++) {
				a.get(y, x) = Ax[i][y];
			}
			cout << "x" << x << ":" << a.getDeterminant()/det <<endl;
		}
	}
	getchar();
}