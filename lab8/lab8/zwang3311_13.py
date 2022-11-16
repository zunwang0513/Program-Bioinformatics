#!/usr/bin/python
#Exercise 13 Reverse String According to Logic Value

string = input("Input your string: ")
val = input("Input logic value (TRUE/FALSE): ")

if val == "TRUE":
    print(string[::-1])
else:
    print(string)
