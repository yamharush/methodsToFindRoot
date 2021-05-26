'''
Authors: yam harush 318886058
         michael harush 208829085
'''
import sympy as sp
from sympy import *
from numpy import *
import numpy as np

#The first method
def derivative(f, x):
    my_f1 = sp.diff(f, x)
    return lambdify(x, my_f1)

def Bisection_Method(func, start, end, epsilon, func_or_derive):
    max_itr = -(np.log(epsilon / (end - start))) / np.log(2)
    iteration = 0
    if func_or_derive == 'fun':
        if func(start) * func(end) >= 0:
            return
    else:
        if func(start) * func(end) > 0:
            print("Bisection method doesn't fit")
            return
    c = start
    while (end- start) > epsilon:
        if iteration > max_itr:
            break
        iteration += 1
        # Find the middle point
        c = (start + end) / 2
        if np.any(np.absolute(c) < epsilon):
            c = 0
        if func(c) == 0.0:
            break
        # Decide the side to repeat the steps
        if func(c) * func(start) < 0:
            end= c
        else:
            start = c
    return c, iteration

def func1_help(func, start_point, end_point):
    iterations = 0
    fx = func
    x1 = sp.symbols('x')
    func = lambdify(x1, func)
    length = 0.1
    s = start_point
    while s < end_point:
        if func(s) * func(s + length) <= 0:
            c, iteration1 = Bisection_Method(func, s, s + length, 0.0001, 'fun')
            iterations += iteration1
            if c is not None:
                print("\n" + "The root is : ", "%.5f" % c)
                print("The number of iterations:", iterations)
                iterations = 0
        s = s + length

    ## derivative check
    s = start_point
    fd = derivative(fx, x1)
    while s < end_point:
        a = fd(s)
        b = fd(s + length)
        if (a * b) <= 0:
            c, iteration2 = Bisection_Method(fd, s, s + length, 10 ** (-10), 'der')
            iterations += iteration2
            if c is not None:
                if func(c) == 0:
                    print("\nThe value of root is : ", "%.5f" % c)
                    print("The number of iterations:", iterations)
                iterations = 0
        s = s + length


# The second method
def Newton_Raphson(f, start, end, epsilon):
    x0 = (end+start)/2
    fx = f
    x1 = start
    x = symbols('x')
    f = lambdify(x, f)
    fd = derivative(fx, x)
    count = 0
    while abs(x1-x0) > epsilon:
        tmp = x1
        if fd(x0) == 0:
            print("division by zero")
            break
        x1 = x0-(f(x0)/fd(x0))
        count += 1
        x0 = tmp
    return x1, count


# The third method
def secant(f, start, end, epsilon):
    x2 = (start+ end) / 2
    x1 = end
    x0 = start
    x = sp.symbols('x')
    f = lambdify(x, f)
    count = 0
    while abs(x2 - x1) > epsilon:
        tmp = x2
        x2 = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0))
        x0 = x1
        x1 = tmp
        count += 1
    return x2, count


def func2_help(f, start, end, message):
    length = 0.1
    x1 = symbols('x')
    fx = lambdify(x1, f)
    x = start
    my_f1 = sp.diff(f, x1)
    f1 = lambdify(x1, my_f1)
    if message == 2:
        while x < end:
            if fx(x) * fx(x + length) <= 0:
                c, count = Newton_Raphson(f, x, x + length, 0.0001)
                if c is not None:
                    print("\nNewton-Raphson num of iterations", count)
                    print("Newton :The value of root is : ", "%.4f" % c, "\n")
            x = x+length
        ## derivative check

        s = start
        while s < end:
            if f1(s) * f1(s + length) <= 0:
                c, count = Newton_Raphson(my_f1, s, s + length, 0.0001)
                if c is not None:
                    if f1(c) == 0:
                        print("Newton-Raphson num of iterations", count)
                        print("Newton :The value of root is : ", "%.4f" % c, "\n")
            s = s + length

    elif message == 3:
        while x < end:
            if fx(x) * fx(x + length) <= 0:
                c, count = secant(f, x, x + length, 0.0001)
                if c is not None:
                    print("\nsecant num of iterations", count)
                    print("Secant :The value of root is : ", "%.4f" % c, "\n")
            x = x + length

        s = start
        while s < end:
            if f1(s) * f1(s + length) <= 0:
                c, count = secant(my_f1, s, s + length, 0.0001)
                if c is not None:
                    if f1(c) == 0:
                        print("secant num of iterations", count)
                        print("Newton :The value of root is : ", "%.4f" % c, "\n")
            s = s + length
    return


def main():
    choice = int((input("Please choose one of the methods to calculate:\n1.Bisection\n2.Newton\n3.Secant\n4.Exit\n")))
    while true:
        x = sp.symbols('x')
        f = x ** 4 + x ** 3 - 3 * x ** 2
        a = -5
        b = 10
        if choice == 4:
            break
        if choice == 1:
            print("\n----------------------Bisection Method--------------------------")
            print("\n~~~~~~~~The first example~~~~~~~~~~~~\n\nF(x)=", f)
            print("The range is: [{},{}]".format(a, b))
            func1_help(f, a, b)
            break
        elif choice == 2:
            print("\n-------------------------Newton-Raphson Method-------------------")
            print("\n~~~~~~~~The first example~~~~~~~~~~~~\n\nF(x)=", f)
            print("The range is: [{},{}]".format(a, b))
            func2_help(f, a, b, 2)
            break
        elif choice == 3:
            print("\n---------------------Secant Method---------------------------")
            print("\n~~~~~~~~The first example~~~~~~~~~~~~\n\nF(x)=", f)
            print("The range is: [{},{}]".format(a, b))
            func2_help(f, a, b, 3)
            break
        else:
            choice = int(input("illegal choice! please try again\n "))


main()
