#!/usr/bin/python
#Exercise 19 mean of columns and mean of mean

mymtrix = [[1,2,3,4,5], [6,7,8,9,10],[11,12,13,14,15]]
print("My matrix is: ")
print(mymtrix)

sumofcolumn = [0] * 5
meanofcolumn = []
meanofmean = 0

for row in mymtrix:
    for columnnum in range(len(row)):
        sumofcolumn[columnnum] += columns[columnnum]
	    
for row2 in range(len(sumofcolumn)):
    meanofcolumn.append(sumofcolumn[row2] / 3)

print("Means of column:")
print(meanofcolumn)

for num in meanofcolumn:
    meanofmean += num
meanofmean /= len(mymtrix[0])
print("Mean of means:")
print(meanofmean)
