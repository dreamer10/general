#include "Solution.h"
#include <iostream>
#include <string>
#include <vector>

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
	    for (int i = 0; i < n; i++) {
            if (A[i] >= 1 && A[i] <= n)
                swap(A[i], A[A[i] - 1]);
        }
		::printArray(A, n);
        
        for (int i = 0; i < n; i++) {
            if (A[i] != i + 1)
			return i + 1;
		}
		return n + 1;
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
	cout << solution.firstMissingPositive(A, 4) << endl;
}
