#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using std::vector;
using std::string;
using std::unordered_map;

struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode(int x) : val(x), left(NULL), right(NULL) {}
};

class Solution {
public:
	Solution() = default;
	int firstMissingPositive(int A[], int n);

	TreeNode *buildTree(vector<int> &preorder, vector<int> &inorder,
										 unordered_map<int, int> &val2Index,
										 int preLo, int preHi,
										 int inLo, int inHi);
	TreeNode *buildTree(vector<int> &preorder, vector<int> &inorder);
	double pow(double x, int n);
	int lengthOfLastWord(const char *s);
	vector<vector<int> > permuteUnique(vector<int> &num);

private:
	void swap(int &a, int &b);
	void permuteUnique(vector<int> &num, vector<int> &tmp, bool canUse[], vector<vector<int> > &result);

};


#endif
