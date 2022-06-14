from statistics import mean


class Student:
    def __init__(self, name: str, courses: "list[str]", marks: "list[float]"):
        self.name = name
        self.courses = courses
        self.marks = marks

    def addCourse(self, course: str, mark: float) -> None:
        """
        Add a course name and the corresponding mark to the appropriate attributes.
        This function doesn't return anything.
        """
        self.course = course
        self.mark = mark

    def meanMark(self) -> float:
        """
        Return the average of a student's marks
        """
        return sum(self.marks)/len(self.marks)

    def transcript(self) -> dict:
        """
        Return a dictionary where the keys are the course names and the values
        are the corresponding marks
        """
        return dict( zip(self.courses, self.marks) )



class Course:
    def __init__(self, name: str, students: "list[Student]" = []):
        self.name = name
        self.students = students

    def addStudent(self, student: Student):
        """
        Check that the student is taking the course and add it to the course list.
        Otherwise raise the error.
        """
        courselist = []
        if self in student.courses:
            courselist.append(student)
        raise RuntimeError(f"{student.name} is not taking {self.name}") # keep


    def meanMark(self) -> float:
        """
        Calculate and return the avarage mark for this course over all enrolled students.
        """
        listofmarks = []
        for student in self.students:
           listofmarks.append(student.transcript()[self.name])
        return mean(listofmarks)


# some tests:  don't modify ----------------------------------------------------
import unittest

anna = Student("Anna Schuppe", ["Maths", "Python", "Science"], [8.5, 7.8, 9.1])
eloy = Student("Eloy Vallina", ["Maths", "Python", "Geography"], [7.8, 9.6, 4.5])
joha = Student("Johanna Zipf", [], [])
py101 = Course("Python", students = [anna, eloy])
sc101 = Course("Science", students = [anna])

class CorrectNumberStudents(unittest.TestCase):
    def test0(self):
        self.assertEqual(len(py101.students), 2)

class StudentNotEnrolled(unittest.TestCase):
    def test0(self):
        with self.assertRaises(RuntimeError):
            sc101.addStudent(eloy)

class StudentMean(unittest.TestCase):
    def test0(self):
        self.assertEqual(eloy.meanMark(), sum([7.8, 9.6, 4.5]) / 3)

class CourseMean(unittest.TestCase):
    def test0(self):
        self.assertEqual(py101.meanMark(), sum([7.8, 9.6]) / 2)

class AddCourse(unittest.TestCase):
    def test0(self):
        joha.addCourse("Physics", 9.5)
        self.assertEqual(len(joha.courses), 1)

if __name__ == '__main__':
    unittest.main()
