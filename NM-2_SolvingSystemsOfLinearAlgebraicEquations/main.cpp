#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#define ANSI_COLOR_BLUE "\033[34m"
#define ANSI_COLOR_RESET "\033[0m"
#define ANSI_COLOR_GREEN "\033[32m"
#define ANSI_COLOR_RED "\033[31m"

using namespace std;

vector<vector<double>> GetMinor(const vector<vector<double>>& matrix, const int& row, const int& column) {

	vector<vector<double>> minor;
	for (int i = 0; i < matrix.size(); i++) {
		if (i == row) continue;
		vector<double> rowMinor;
		for (int j = 0; j < matrix[i].size(); j++) {
			if (j == column) continue;
			rowMinor.push_back(matrix[i][j]);
		}
		minor.push_back(rowMinor);
	}
	return minor;
}

double GetDeterminate(const vector<vector<double>>& matrix) {

	if (matrix.size() == 1) return matrix[0][0];
	if (matrix.size() == 2) {
		double a = matrix[0][0]; double b = matrix[0][1];
		double c = matrix[1][0]; double d = matrix[1][1];
		return a * d - c * b;
	}
	double det = 0;
	for (int j = 0; j < matrix.size(); j++) {
		double coef = matrix[0][j];
		double sign;
		if (j % 2 == 0) sign = 1;
		else sign = -1;
		det += sign * coef * GetDeterminate(GetMinor(matrix, 0, j));
	}
	return det;
}

vector<vector<double>> ReplaceMatrixColumn(const vector<vector<double>>& matrix, const vector<double>& B, const int& column) {

	vector<vector<double>> newMatrix(matrix.size(), vector<double>(matrix[0].size(), 0));
	int index = 0;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			if (j == column)  { newMatrix[i][j] = B[index]; index++; }
			else newMatrix[i][j] = matrix[i][j];
		}
	}
	return newMatrix;
}

void CramerMethod(const vector<vector<double>>& matrix, const vector<double>& B) {

	if (matrix.empty() || matrix[0].size() != matrix.size() || B.size() != matrix.size()) {
		throw "Error: <Invalid matrix size>";
	}
	double det = 0;
	vector<double> detX;
	vector<double> x;
	det = GetDeterminate(matrix);
	if (abs(det) < 1e-9) throw "Error: <System has no unique solution (det(A) = 0)";
	for (int i = 0; i < matrix.size(); i++) {
		detX.push_back(GetDeterminate(ReplaceMatrixColumn(matrix, B, i)));
	}
	cout << "detXn : ";
	for (auto i : detX) {
		cout << i << " ";
	}cout << endl;
	for (int i = 0; i < matrix.size(); i++) {
		x.push_back(detX[i] / det);
	}
	cout << "Xn : ";
	for (auto i : x) {
		cout << i << " ";
	}cout << endl;
}

void MatrixShow(const vector<vector<double>>& matrix) {

	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[0].size(); j++) {
			cout << matrix[i][j] << " ";
		}cout << endl;
	}
}

vector<vector<double>> GetInvertedMatrix(const vector<vector<double>>& matrix) {

	double sign;
	vector<vector<double>> matrixComplement(matrix.size(), vector<double>(matrix[0].size(), 0));
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[0].size(); j++) {
			if ((i + j) % 2 == 0) sign = 1;
			else sign = -1;
			matrixComplement[i][j] = sign * (GetDeterminate(GetMinor(matrix, i, j)));		}
	}
	cout << "Matrix : \n";
	MatrixShow(matrixComplement);
	for (int i = 0; i < matrixComplement.size(); i++) {
		for (int j = i + 1; j < matrixComplement[0].size(); j++) {
			swap(matrixComplement[i][j], matrixComplement[j][i]);
		}
	}
	cout << "Matrix : \n";
	MatrixShow(matrixComplement);
	double det = GetDeterminate(matrix);
	for (int i = 0; i < matrixComplement.size(); i++) {
		for (int j = 0; j < matrixComplement[0].size(); j++) {
			matrixComplement[i][j] = matrixComplement[i][j] / det;
		}
	}
	cout << "Matrix : \n";
	MatrixShow(matrixComplement);

	return matrixComplement;
}

void MatrixMethod(const vector<vector<double>>& matrix, const vector<double>& B) {

	if (matrix.empty() || matrix[0].size() != matrix.size() || B.size() != matrix.size()) {
		throw "Error: <Invalid matrix size>";
	}
	vector<vector<double>> tmpMatrix;
	vector<double> Xn(matrix.size(), 0);
	double det = GetDeterminate(matrix);
	if (abs(det) < 1e-9) throw "Error: <System has no unique solution (det(A) = 0)>";
	tmpMatrix = GetInvertedMatrix(matrix);
	for (int i = 0; i < tmpMatrix.size(); i++) {
		int bIndex = 0;
		for (int j = 0; j < tmpMatrix[0].size(); j++) {
			Xn[i] += tmpMatrix[i][j] * B[bIndex];
			bIndex++;
		}
	}
	cout << "Xn : ";
	for (auto i : Xn) {
		cout << i << " ";
	}cout << endl;
}

void GaussForward(vector<vector<double>>& matrix) {

	for (int i = 0; i < matrix.size(); i++) {
		int maxRow = i;
		for (int j = i + 1; j < matrix.size(); j++) {
			if (abs(matrix[j][i] > abs(matrix[maxRow][i]))) {
				maxRow = j;
			}
		}
		swap(matrix[i], matrix[maxRow]);
		double fixPlace = matrix[i][i];
		if (abs(fixPlace) < 1e-9) throw "Error: <System has no unique solution (det(A) = 0)>";
		for (int j = i; j <= matrix.size(); j++) {
			matrix[i][j] /= fixPlace;
		}
		for (int m = i + 1; m < matrix.size(); m++) {
			double fixCell = matrix[m][i];
			for (int j = i; j <= matrix.size(); j++) {
				matrix[m][j] -= matrix[i][j] * fixCell;
			}
		}
	}	
}

vector<double> GaussReverse(const vector<vector<double>>& matrix) {

	vector<double> Xn(matrix.size(), 0);
	for (int i = matrix.size() - 1; i >= 0; i--) {
		Xn[i] = matrix[i][matrix.size()];
		for (int j = i + 1; j < matrix.size(); j++) {
			Xn[i] -= matrix[i][j] * Xn[j];
		}
	}
	return Xn;
}

void GaussMethod(const vector<vector<double>>& matrix, const vector<double>& B) {

	if (matrix.empty() || matrix[0].size() != matrix.size() || B.size() != matrix.size()) {
		throw "Error: <Invalid matrix size>";
	}
	vector<vector<double>> matrixExpanded(matrix.size(), vector<double>(matrix[0].size() + 1, 0));
	vector<double> Xn;
	int index = 0;
	cout << matrix.size() << endl;
	for (int i = 0; i < matrixExpanded.size(); i++) {
		for (int j = 0; j < matrixExpanded[0].size(); j++) {
			if (j == matrixExpanded[0].size() - 1) { matrixExpanded[i][j] = B[index]; index++; }
			else matrixExpanded[i][j] = matrix[i][j];
		}
	}
	GaussForward(matrixExpanded);
	MatrixShow(matrixExpanded);
	Xn = GaussReverse(matrixExpanded);
	cout << "Xn : ";
	for (auto i : Xn) {
		cout << i << " ";
	}cout << endl;
}

int main() {

	vector<vector<double>> matrix2 = {
		{ 1, -1, -1},
		{-2, 1, 2},
		{3, 3, -1}
	};
	vector<double> B2 = { 0, -1, 8 };

	vector<vector<double>> matrixTest = {
		{0.79, 0.05, -0.25, 0.08},
		{0.21, -0.13, 0.27, -0.8},
		{-0.11, -0.84, 0.35, 0.06},
		{-0.08, 0.15, -0.5, -0.12}
	};
	vector<double> B = { 2.15, 0.44, -0.83, 1.16 };

	try {
		//MatrixMethod(matrixTest, B);
	}
	catch (char* err) {
		cerr << err << endl;
	}
	try {
		//CramerMethod(matrixTest, B);
	}
	catch(char* err) {
		cerr << err << endl;
	}	
	try {
		GaussMethod(matrixTest, B);
	}
	catch (char* err) {
		cerr << err << endl;
	}

	return 0;
}