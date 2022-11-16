#!/usr/bin/python
#Exercise 16 sum and product of two nums in a list

rawstr = input("Give a list of numbers separated by space: ")
lst = [int(i) for i in rawstr.split()]

rawposition = input("Give two integers of position separated by space (First number position is 1): ")
position = [int(j) for j in rawposition.split()]

sumoftwo = lst[position[0]-1] + lst[position[1]-1]
productoftwo = lst[position[0]-1] * lst[position[1]-1]
print("sum is " + str(sumoftwo))
print("product is " + str(productoftwo))
