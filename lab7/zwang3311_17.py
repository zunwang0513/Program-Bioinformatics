#!/usr/bin/python
#Exercise 17 Find median of two lists

flst = input("Give me list one separated by space: ")
slst = input("Give me list two separated by space: ")

lst1 = [int(i) for i in flst.split()]
lst2 = [int(i) for i in slst.split()]

sortedlst1 = sorted(lst1)
sortedlst2 = sorted(lst2)

idx1 = (len(lst1) - 1) // 2
idx2 = (len(lst2) - 1) // 2

output = []
output.append(sortedlst1[idx1])
output.append(sortedlst2[idx2])

print(output)
