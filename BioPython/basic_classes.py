#!/usr/bin/env python
#basic_classes.py
import sys

class circle():
#defining class with some attributes
    def __init__(self,color,radius):
        self.color=str(color)
        self.radius=float(radius)
#defining the method
    def diameter(self):
        return 2*(self.radius)
    def circumference(self):
        return 2*3.14*(self.radius)
    def isRed(self):
        if self.color=="red":
          return True
class GraduateStudent():
    def __init__(self,first_name,last_name,year,major):
        self.first_name=first_name
	self.last_name=last_name
	self.year=year
        self.major=major
    def year_matriculated(self):
        return 2020-self.year
if __name__ == "__main__":
    Circle1 = circle("red",2)
    Circle2 = circle("Blue",4)
    print("the Diameter, Circumference and 'IS color Red?' for the first circle is ", Circle1.diameter(), Circle1.circumference(), Circle1.isRed())

    print("the Diameter, Circumference and 'IS color Red?' for the second circle is ", Circle2.diameter(), Circle2.circumference(), Circle2.isRed())
        
    year2015 = GraduateStudent("dsg", "reg", 1, "ergs")

    print("The year matriculated is", year2015.year_matriculated())    
      

