#!/usr/bin/env python3
#1. multiply.py

def multiply(a,b):
    """
    Returning the multiplication values of two given numbers a and b
    """
    return a*b


#2 hello_name
def hello_name(my_name="you"):
    """
    settting the degault name as you
    """
    return("Hello, {}!".format(my_name))

#3 less_than_ten

def less_than_ten(list1):
    """
    Returns the least of the 10 digits for a list of numbers
    """
    new_list=[]

    for x in list1:
        if x<10:
            new_list.append(x)
    return(new_list)

