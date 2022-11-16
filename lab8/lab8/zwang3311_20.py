#!/usr/bin/python
#Exercise 20 pairwise distance

# I don't know if I understand the problem correctly.
# I computed the sum of pairwise difference of each element compared
# to all other elements in the matrix and calculate its average
# Resulting in a same size 3 x 5 matrix

my_mtrix = [[1,2,3,4,5], [6,7,8,9,10], [11,12,13,14,15]]
print("This is my matrix:")
print(my_mtrix)

diff = [[0] * 5] * 3
avgdiff = [[0] * 5] * 3
for m in range(len(my_mtrix)):
    for n in range(len(my_mtrix[0])):
        for i in range(len(my_mtrix)):
            for j in range(len(my_mtrix[0])):
                diff[m][n] = diff[m][n] + my_mtrix[m][n] - my_mtrix[i][j]
print(diff)

for i in range(len(diff)):
    print(i)
    for j in range(len(diff[0])):
        print(j)
        print(diff[i][j])
        avgdiff[i][j] = diff[i][j] / 14
        print(avgdiff[i][j])
    #Since there are 14 pairwise differences, did not count difference with self

print(avgdiff)
