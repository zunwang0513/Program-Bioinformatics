#!/usr/bin/python
#Exercise 20 pairwise distance

# I don't know if I understand the problem correctly.
# I computed the sum of pairwise difference of each element compared
# to its right neighbour in the row and calculate average of all number in the matrix
# Resulting in a same size 3 x 4 pairwise matrix and an interger average

my_mtrix = [[1,2,3,4,5], [6,7,8,9,10], [11,12,13,14,15]]
print("This is my matrix:")
print(my_mtrix)

diff = [[0] * 4] * 3

for i in range(len(my_mtrix)):
	for j in range(len(my_mtrix[0]) - 1):
		diff[i][j] += my_mtrix[i][j] - my_mtrix[i][j + 1]
sumofnum = 0
for m in range(len(diff)):
	for n in range(len(diff[0])):
		sumofnum += diff[i][j]
average = sumofnum / 12
print("Pairwise average is ", average)
