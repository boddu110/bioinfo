#!/usr/bin/env python3
import sys

def hello_name (str="you!"):
    

    return("Hello, {}".format(str))

if __name__ == "__main__":
    arg_count = len(sys.argv) - 1
    if arg_count != 1:
        print(hello_name())
    else:
        print("Hello, {}!".format(sys.argv[1]))

