#!/usr/bin/python
#Exercise 9 Return sum unless equal

a = int(input("Input a: "))
b = int(input("Input b: "))

if a == b:
    print (a + b + a * a + b * b)
else:
    print (a + b)
