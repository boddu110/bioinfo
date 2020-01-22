#!/usr/bin/env python
#fibbonacci.py
import sys
def population(n,k):
    if n ==0:
         return 0
    elif n ==1:
        return 1

    if n >1000:
        sys.exit("The number  need to be between 1 and 10000 ")
    elif  k>10000:
        sys.exit("The number should be between 1 and 10000 ")

    population = [0] * n
    population[0] = 1
    population[1] = 1

    for i in range(2, n):
        offspring = population[i-2] * k
        population[i] = population[i-1] + offspring

    return population[n-1]

if __name__ == '__main__':
     arg_count = len(sys.argv) -1
     if arg_count < 2:
         raise Exception("This Script required 2 arguments")
     a1,a2 = [int(x) for x in sys.argv[1:3]]
     pop = population(a1,a2)

     print(pop)
                  
