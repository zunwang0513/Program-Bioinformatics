#!/usr/bin/python
#Excercise 7 Swapping values

a = 10
b = 30

print("a is 10 and b is 30. Now Swap 1 using temp")
temp = a
a = b
b = temp
print ("a is now " + str(a))
print ("b is now " + str(b))

a = 10
b = 30
print("a is 10 and b is 30. Now Swap 2 without temp")
a, b = b, a
print ("a is now " + str(a))
print ("b is now " + str(b))
