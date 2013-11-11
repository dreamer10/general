#include "Solution.h"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

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

int main(int argc, char *argv[])
{
	int A[] = {3, 4, -1, 1};
	Solution solution;
	vector<int> a = {-1};
	vector<int> b = {-1};
	TreeNode *root = solution.buildTree(a, b);
	cout << root->val << endl;
}
