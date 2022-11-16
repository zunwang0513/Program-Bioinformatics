#!/usr/bin/python
#Excercise 6 Left Pascal Triangle

x = 1
y = 9
star = x
while y > 0:
    if x <= 5:
        star = x
    else:
        star = y
    while star > 0:
        print("*",end = '')
        star -= 1
    print()
    x += 1
    y -= 1s
