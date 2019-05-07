# test.jl

Most functions in **PosDefManifold** are tested, both for real and complex data input. This unit declares the function `testall()` that performs all tests.

Some functions are fully tested, the others are just executed.
Unce you ran it, for each method of each function,
a ⭐ sign is printed if the test is succesful, while
a ⛔ sign is printed if the test is not succesful.
A ☆ sign is printed if the function has been executed correctly.

If there are fails, the concerned functions will be listed as Warnings.

The first time you execute the test it will take some time.
