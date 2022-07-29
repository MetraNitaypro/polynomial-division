#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
uint64_t** divM(uint64_t *a, uint64_t *b, int size);//деление многочленов 
void deg(uint64_t* c, int size, int* deg);//определение степени многочленов
void shift(uint64_t* res, uint64_t* input, int size, int shift_1, int shift_2, bool tag, int subsum);//сдвиг
void print(uint64_t* a, int size, int l);
void koder(uint64_t* res, uint64_t* m, uint64_t* g, int size, int r);
bool dekoder(uint64_t* a, uint64_t* g, int size);
vector<vector<uint64_t>> GenVec(int k, int l, int d_1, int deg_g, uint64_t g, ofstream & f_1);
int c(int m, int n);
int d(vector<uint64_t>& codeWord, int size);
int w(uint64_t val, int n);
void Gen_e(vector<uint64_t>& e, uint64_t currentVal, int m, int n, int current_count, int countBit);
int main() {
	//ifstream f_1("input.txt");
	//int k = 0;
	//int l = 0;
	//int n = 0;
	//f_1 >> k;
	//f_1 >> l;
	//n = (l + k - 1) / (sizeof(uint64_t) * 8);
	//if ((l + k - 1) % (sizeof(uint64_t) * 8) != 0) {
	//	n++;
	//}
	//uint64_t* m = new uint64_t[n];
	//uint64_t* g = new uint64_t[n];
	//uint64_t* a = new uint64_t[n];
	//uint64_t* e = new uint64_t[n];
	//memset(m, 0, sizeof(uint64_t) * n);
	//memset(g, 0, sizeof(uint64_t) * n);
	//memset(a, 0, sizeof(uint64_t) * n);
	//memset(e, 0, sizeof(uint64_t) * n);
	//int count = k % (sizeof(uint64_t) * 8);
	uint64_t val = 0;
	//int i = k / (sizeof(uint64_t) * 8);

	//for (; i >= 0 ; i--) {
	//	for (int j = count - 1; j >= 0 ; j--) {
	//		f_1 >> val; 
	//		g[i] |= val << j;
	//	}
	//	count = (sizeof(uint64_t) * 8) - 1;
	//}

	//i = l / (sizeof(uint64_t) * 8);
	//count = l % (sizeof(uint64_t) * 8);
	//for (; i >= 0; i--) {
	//	for (int j = count - 1; j >= 0; j--) {
	//		f_1 >> val;
	//		m[i] |= val << j;
	//	}
	//	count = (sizeof(uint64_t) * 8) - 1;
	//}

	//i = n / (sizeof(uint64_t) * 8);
	//count = (l + k - 1) % (sizeof(uint64_t) * 8);
	//for (; i >= 0; i--) {
	//	for (int j = count - 1; j >= 0; j--) {
	//		f_1 >> val;
	//		e[i] |= val << j;
	//	}
	//	count = (sizeof(uint64_t) * 8) - 1;
	//}
	//cout << "g(x) ";
	//print(g, n);
	//cout << "m(x) ";
	//print(m, n);
	//
	//koder(a, m, g, n, k - 1);
	////print(a, n);
	//cout << "e(x) ";
	//print(e, n);
	//for (int kk = 0; kk < n; kk++) {
	//	a[kk] ^= e[kk];
	//}
	//cout << "b(x)" << endl;
	//print(a, n);
	//dekoder(a, g, n);
	uint64_t G = 0;
	int count = 4;
	count = 4;
	int SizeG;
	cin >> SizeG;
	for (int j = SizeG - 1; j >= 0; j--) {
		cin >> val;
		G |= val << j;
	}
	vector<vector<vector<uint64_t>>> result;
	ofstream f_1("results.txt");
	for (int ii = 0; ii < 3; ii++) {
		for (int j_1 = 3; j_1 < 4; j_1++) {
			 GenVec(4, ii, j_1, SizeG, G, f_1);
		}
	}

	f_1.close();
	cin >> SizeG;
	for (int j = SizeG - 1; j >= 0; j--) {
		cin >> val;
		G |= val << j;
	}
	ofstream f_2("results_1.txt");
	for (int ii = 0; ii < 3; ii++) {
		for (int j_1 = 4; j_1 < 5; j_1++) {
			GenVec(3, ii, j_1, SizeG, G, f_2);
		}
	}
	
	f_2.close();

	return 0;
}

void Gen_e(vector<uint64_t>& e, uint64_t currentVal, int m, int n,  int current_count, int countBit) {
	if (countBit == m) {
		e.push_back(currentVal);
		return;
	}
	for (int i = current_count; i < n; i++) {
		currentVal |= 1 << i;
		Gen_e(e, currentVal, m, n, ++current_count, ++countBit);
		currentVal ^= 1 << i;
		countBit--;
	}
	return;
}
int c(int m, int n) {
	int k = 1;
	for (int i = 0; i < n; i++) {
		k *= i + 1;
	}
	for (int i = 0; i < m; i++) {
		k /= i + 1;
	}
	for (int i = 0; i < (m - n); i++) {
		k /= i + 1;
	}
	return k;
}
int w(uint64_t val, int n) {
	int k = 0;
	uint64_t mask = 0;
	for (int i = n - 1; i >= 0; i--) {
		mask = 1 << i;
		if (val & mask != 0) {
			k++;
		}
	}
	return k;
}
int d(vector<uint64_t> & codeWord, int size) {
	int d = 10000;
	int count_1 = 0;
	uint64_t mask = 0;
	uint64_t val = 1;
	for (int i = 0; i < codeWord.size() - 1; i++) {
		for (int j = i + 1; j < codeWord.size(); j++) {
			for (int k = 0; k < sizeof(uint64_t) * 8; k++) {
				mask =  val << k;
				if ((codeWord[i] & mask) != (codeWord[j] & mask)) {
					count_1++;
				}
				mask = 0;
				val = 1;
			}
			if (count_1 < d) {
				d = count_1;
				//cout << "d  " << d << "   " << i <<"    " << j << endl;
			}
			count_1 = 0;
		}
	}
	return d;
}
vector<vector<uint64_t>> GenVec(int k, int l, int d_1, int deg_g, uint64_t g, ofstream& f_1) {
	cout << "k = " << k << "  l = " << l << "  d = " << d_1 << "  deg_g = " << deg_g - 1 << "  g = ";
	print(&g, 1, deg_g - 1);
	int size = 1;
	vector<uint64_t> res(1);
	vector<uint64_t> codeW(1);
	vector<uint64_t> codeW_1(1);
	vector<uint64_t> e;
	vector<uint64_t> result;
	int* deg_1 = new int[2];
	int* deg_2 = new int[2];
	deg(&g, size, deg_1);
	for (int j = 1; j < 2; j++) {
		res[0] = j;
		deg(&res[0], size, deg_2);
		koder(&codeW[0], &res[0], &g, size, deg_1[1]);
	}
	int d_2 = 0;
	d_2 = d(codeW, size);
	vector<vector<uint64_t>> res_1(3);
	for (int j = 1; j < d_1; j++) {
		Gen_e(e, 0, j, l + k + deg_1[1], 0, 0);	
	}
	cout << "a " << endl;
	print(&codeW[0], 1, l + k + deg_1[1] - 1);
	int count_10 = 0;
	for (int j = 0; j < 1; j++) {
		for (int jj = 0; jj < e.size(); jj++) {
			codeW_1[j] = codeW[j] ^ e[jj];
			if (dekoder(&codeW_1[j], &g, size)) {
				cout << "e " << endl;
				print(&e[jj], 1, l + k + deg_1[1] - 1);
				cout << "b " << endl;
				print(&codeW_1[j], 1, l + k + deg_1[1] - 1);
				count_10++;
			}
		}
	}
	f_1 << l << " " << count_10 << endl;
	res_1[2].push_back(d_1);
	delete[] deg_1;
	delete[] deg_2;
	return res_1;
}
void print(uint64_t* a, int size, int l) {
	int* deg_1 = new int[2];
	deg(a, size, deg_1);
	if (deg_1[0] == -1) {
		cout << 0 << endl;
		return;
	}
	int count = l;
	uint64_t mask = 1;
	for (int j = deg_1[0]; j >= 0; j--) {
		for (int h = count; h >= 0; h--) cout << ((a[j] >> h) & mask) << " ";
		count = sizeof(uint64_t) * 8 - 1;
	}
	cout << endl;
	delete[] deg_1;
}
uint64_t** divM(uint64_t* a, uint64_t* b, int size) {// деление многочленов
	//cout << "division" << endl;
	uint64_t* a_1 = new uint64_t[size];
	uint64_t* b_1 = new uint64_t[size];
	memset(b_1, 0, sizeof(uint64_t) * size);
	memcpy(a_1, a, sizeof(uint64_t) * size);
	int* deg_1 = new int [2];
	int* deg_2 = new int[2];
    deg(a, size, deg_1);
    deg(b, size, deg_2);
	int subsum = deg_1[0] * sizeof(uint64_t) * 8 + deg_1[1] - (deg_2[0] * sizeof(uint64_t) * 8 + deg_2[1]);
	uint64_t** arr = new uint64_t * [2];
	arr[0] = new uint64_t[size];//частное
	arr[1] = new uint64_t[size];//остаток
	memset(arr[0], 0, sizeof(uint64_t) * size);
	memset(arr[1], 0, sizeof(uint64_t) * size);
	if (deg_1[1] == -1 || deg_2[1] == -1) {
		if (deg_2[1] == -1) {
			delete[] arr;
			return NULL;
		}
		memcpy(arr[1], a_1, sizeof(uint64_t) * size);
		return arr;
	} 
	else if (subsum < 0) {
		memcpy(arr[1], a_1, sizeof(uint64_t) * size);
		return arr;
	}
	else {
		while (subsum >= 0) {
			arr[0][deg_1[0] - deg_2[0]] |= (uint64_t) 1 << (deg_1[1] - deg_2[1]);
			shift(a_1, b, size, deg_1[1] - deg_2[1], deg_1[0] - deg_2[0], true, subsum);
			deg(a_1, size, deg_1);
			subsum = deg_1[0] * sizeof(uint64_t) * 8 + deg_1[1] - (deg_2[0] * sizeof(uint64_t) * 8 + deg_2[1]);
			//print(a_1, size);
		}
		memcpy(arr[1], a_1, sizeof(uint64_t) * size);
	}
	delete[] deg_1;
	delete[] deg_2;
	delete[] a_1;
	delete[] b_1;
	return arr;
}
void deg(uint64_t* c, int size, int* deg) {// определение степени многочлена
    int count = sizeof(uint64_t) * 8 - 1;
	int i = 0;
	for (i = size - 1; i >= 0; i--) {
		while (c[i] >> count == 0) {
			count--;
			if (count < 0)  break;
		}
		if (count >= 0)  break;
		count = sizeof(uint64_t) * 8 - 1;
	}
	deg[0] = i;
	deg[1] = count;
}
void shift(uint64_t* res, uint64_t* input, int size, int shift_1, int shift_2, bool tag, int subsum) {
	if ((subsum % (sizeof(uint64_t) * 8 * 2)) < (sizeof(uint64_t) * 8) && (shift_2 > 0)) shift_2--;
	shift_1 = shift_1 % (sizeof(uint64_t) * 8);
	uint64_t* inputcopy = new uint64_t[size];
	memset(inputcopy, 0, sizeof(uint64_t) * size);
	memcpy(&(inputcopy[shift_2]), input, sizeof(uint64_t) * (size - shift_2));
	if(!tag) memset(res, 0, sizeof(uint64_t) * size);
	for (int i = 0; i < size; i++) {
		if (i == 0) {
			uint64_t h = inputcopy[i] << shift_1;
			res[i] ^= h;
		}
		else {
			res[i] ^= ((inputcopy[i] << shift_1) | (inputcopy[i - 1] >> (sizeof(uint64_t) * 8 - shift_1)));
	
		}
	}	
	delete[] inputcopy;
}
void koder(uint64_t* res, uint64_t* m, uint64_t* g, int size, int r) {
	int subsum = r;
	shift(res, m, size, r, r / (sizeof(uint64_t) * 8), false, subsum);
	uint64_t** res_1;
	res_1 = divM(res, g, size);
	for (int i = 0; i < size; i++) {
		res[i] |= res_1[1][i];
	}
	//cout << "a(x)" << endl;
	//print(res, size);
	delete[] res_1;
}
bool dekoder(uint64_t* a, uint64_t* g, int size) {
	uint64_t** res_1;
	res_1 = divM(a, g, size);
	for (int i = 0; i < size; i++) {
		if (res_1[1][i] != 0) {
			//cout << "E = 1 " << "s(x) = ";
			//print(res_1[1], size);
			delete[] res_1;
			return false;
		}
	}
	//cout << "E = 0 " << "s(x) = ";
	//print(res_1[1], size);
	delete[] res_1;
	return true;
}