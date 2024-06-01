#!/usr/bin/env python
# coding: utf-8

# In[6]:



# Algorithm Description and Time Complexity Analysis

# Algorithm Description:
# The algorithm solves the Minimum-Value Grouping Problem using dynamic programming. 
# Given a sequence of n positive integers, it finds the optimal way to partition the sequence into k groups 
# to minimize the sum of the cubes of the sums of the numbers in each group.
# The dynamic programming state dp[i][j] represents the minimum objective function value for partitioning the 
# first i elements into j groups. To find the optimal solution, we consider every possible last group formed by 
# elements x to i, where x is in the range [j-1, i], and find the one that minimizes the objective function.
# The algorithm uses a bottom-up approach to fill in the values of the dp table, iterating over the number of 
# elements and the number of groups to compute the values based on previous computations.

# Time Complexity Analysis:
# The algorithm consists of:
# 1. Calculating prefix sums, which takes O(n) time.
# 2. Initializing DP and cut position tables, which takes O(nk) time.
# 3. Filling the DP table, which takes O(n^2k) time due to the three nested loops.
# 4. Retrieving the cut positions, which takes O(k) time.
# Therefore, the overall time complexity of the algorithm is O(n^2k).

# Note: Since we are given that k ≤ n, in the worst case where k is close to n, the complexity remains bounded by O(n^3), 
# which is a polynomial-time complexity relative to the input size. However, the algorithm may not be efficient for large values 
# of n and k due to the cubic factor.



# This function calculates the minimum value of grouping the input sequence into k groups, 
# such that the sum of the cube of each group's sum is minimized. It uses a dynamic programming approach.
def min_sum_of_cubes(sequence, k):
    n = len(sequence)
    # Compute the prefix sums for the sequence to enable O(1) sum queries.
    # This preprocessing step takes O(n) time.
    prefix_sums = [0] * (n+1)
    for i in range(1, n+1):
        prefix_sums[i] = prefix_sums[i-1] + sequence[i-1]

    # A helper function to calculate the sum of cubes of a segment in the sequence
    # using the prefix sums, which takes constant time, O(1).
    def sum_of_cubes(i, j):
        return (prefix_sums[j] - prefix_sums[i]) ** 3

    # Initialize the dynamic programming table with infinity, which represents
    # the cost for each (i, j) where i is the number of elements and j is the number of groups.
    # The table has (n+1) rows and (k+1) columns, so this initialization takes O(nk) time.
    dp = [[float('inf')] * (k+1) for _ in range(n+1)]
    dp[0][0] = 0  # Base case: no cost for 0 elements in 0 groups

    # Initialize the cut positions table to trace the group divisions.
    # This also takes O(nk) time.
    cut_positions = [[0] * (k+1) for _ in range(n+1)]

    # Fill the dynamic programming table
    # The outer two loops run in O(nk) time since they iterate through each cell in the DP table.
    for i in range(1, n+1):
        for j in range(1, min(i, k)+1):
            # The inner loop tries all possible last cut positions,
            # which takes O(n) time in the worst case for each cell.
            for x in range(j-1, i):
                value = dp[x][j-1] + sum_of_cubes(x, i)
                if value < dp[i][j]:
                    dp[i][j] = value
                    cut_positions[i][j] = x

    # Retrieve the cut positions from the cut_positions table,
    # which takes O(k) time since we have at most k cuts to trace back.
    cuts = [n]
    current = n
    for j in range(k, 0, -1):
        cuts.append(cut_positions[current][j])
        current = cut_positions[current][j]
    cuts = cuts[::-1]

    return dp[n][k], cuts

# Read the input
n = int(input())
k = int(input())
sequence = [int(i) for i in input().split()]

# Solve the problem
optimal_value, cuts = min_sum_of_cubes(sequence, k)

# Print the results
print(optimal_value)
print(cuts)

