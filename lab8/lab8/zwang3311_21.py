#!/usr/bin/python
#Exercise 20 Pairwise Difference

mymtrix = [[1,2,3,4,5], [6,7,8,9,10],[11,12,13,14,15]]
print("Original Matrix: ")
print(mymtrix)

for row in mymtrix:
	row.reverse()
print("Altered Matrix: ")
print(mymtrix)
