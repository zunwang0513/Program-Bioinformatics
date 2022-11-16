#!/usr/bin/python
#Exercise 18 Create Dictionary

num = int(input("Give an integer: "))
dicts = {}
for i in range(1,num + 1):
	dicts[i] =  i * i

print(dicts)
