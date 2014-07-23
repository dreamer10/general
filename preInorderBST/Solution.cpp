#include "Solution.h"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>

using namespace std;

void printArray(int A[], int n)
{
	for (int i = 0; i < n; i++) {
		cout << A[i] << " ";
	}
	cout << endl;
}

int Solution::firstMissingPositive(int A[], int n)
{
	    for (int i = 0; i < n;) {
            if (A[i] >= 1 && A[i] <= n && A[i] - 1 != i && A[i] != A[A[i] - 1])
                swap(A[i], A[A[i] - 1]);
			else
				i++;
        }
		::printArray(A, n);
        
        for (int i = 0; i < n; i++) {
            if (A[i] != i + 1)
			return i + 1;
		}
		return n + 1;
}

TreeNode *Solution::buildTree(vector<int> &preorder, vector<int> &inorder)
{
	unordered_map<int, int> val2Index;

	if (preorder.empty() || inorder.empty())
		return NULL;

	for (int i = 0; i < inorder.size(); i++) {
		int val = inorder[i];
		int iter = i;
		val2Index.insert(std::make_pair<int, int>(std::move(val), std::move(iter)));
	}

	return buildTree(preorder, inorder, val2Index, 0, int(preorder.size()) - 1, 0, int(preorder.size()) - 1);
}

TreeNode *Solution::buildTree(vector<int> &preorder, vector<int> &inorder,
					 unordered_map<int, int> &val2Index, 
					 int preLo, int preHi, int inLo, int inHi)
{
	if (preLo > preHi || inLo > inHi)
		return NULL;

	int rootVal = preorder[preLo];
	TreeNode *node = new TreeNode(rootVal);
	int pivot = val2Index[rootVal];
	cout << "pivot: " << pivot << endl;

	cout << "1" << endl;

	int leftLen = pivot - inLo;
	int rightLen = inHi - pivot;
	
	TreeNode *left = buildTree(preorder, inorder, val2Index, preLo + 1, preLo + leftLen, inLo, pivot - 1);

	cout << "2" << endl;

	TreeNode *right = buildTree(preorder, inorder, val2Index, preLo + leftLen + 1, preHi, pivot + 1, inHi);

	cout << "3" << endl;

	node->left = left;
	node->right = right;
	return node;
}


void Solution::swap(int &a, int &b)
{
	int tmp = a;
	a = b;
	b = tmp;
}

double Solution::pow(double x, int n)
{
	if (n == 0)
		return 1;
	else if (n == 1)
		return x;
	
	if (n < 0) {
		x = 1 / x;
		n = -n;
	}
	if ((n & 1) == 0) {
		double y = pow(x, n >> 2);
		return y * y;
	} else if ((n & 1) == 1) {
		cout << "n: " << n << endl;
		cout << ((n - 1) >> 1) << endl;
		double y = pow(x, (n - 1) >> 1);
		cout << "true" << endl;
		return y * y * x;
	}
}

int Solution::lengthOfLastWord(const char *s) {
	// IMPORTANT: Please reset any member data you declared, as
	// the same Solution instance will be reused for each test case.
	const char *iter = s;
	string str;
	cout << " 1 " << endl;
	while (*iter != '\0') {
	cout << " 2 " << endl;
		str += *iter;
		iter++;
	}
	cout << " 3 " << endl;
	stringstream ss(str);
	cout << " 4 " << endl;
	vector<string> vect;
	cout << " 5 " << endl;
	string line;
	cout << " 6 " << endl;
	while (getline(ss, line, ' ')) {
		cout << " 7 " << endl;
		vect.push_back(line);
	}
	cout << " 8 " << endl;
	if (vect.empty())
			return 0;
	if (vect[vect.size() - 1].empty()) {
		cout << " 9 " << endl;
		return 0;
	} else {
		return vect[vect.size() - 1].size();
	}
}

vector<vector<int> > Solution::permuteUnique(vector<int> &num) {
	// IMPORTANT: Please reset any member data you declared, as
	// the same Solution instance will be reused for each test case.
	std::sort(num.begin(), num.end());
	bool *canUse = new bool[num.size()];
	memset(canUse, true, sizeof(bool));
	vector<vector<int> > result;
	vector<int> tmp;
	permuteUnique(num, tmp, canUse, result);
	delete[] canUse;
	return result;
}

void Solution::permuteUnique(vector<int> &num, vector<int> &tmp, bool canUse[], vector<vector<int> > &result)
{
	if (tmp.size() == num.size()) {
		result.push_back(tmp);
		return;
	}
	for (int i = 0; i < num.size(); i++) {
		if (canUse[i]) {
			if (i != 0 && num[i - 1] == num[i] && canUse[i - 1])
				continue;
			tmp.push_back(num[i]);
			canUse[i] = false;
			permuteUnique(num, tmp, canUse, result);
			tmp.pop_back();
			canUse[i] = true;
		}
	}
}


int main(int argc, char *argv[])
{
	int A[] = {3, 4, -1, 1};
	Solution solution;
	vector<int> a = {-1};
	vector<int> b = {-1};
	//TreeNode *root = solution.buildTree(a, b);
	//cout << root->val << endl;
	//cout << - 2 ^ 32 << endl;
	double x = 8.88023;
	int n = 3;
	char *s = "";
	cout << solution.lengthOfLastWord(s) << endl;

	vector<int> num = {-1, 2, -1, 2, 1, -1, 2, 1};

}
