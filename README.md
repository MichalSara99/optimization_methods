# Optimization Methods
This is a header-only C++ library containing several well-known optimization methods. 
The library is written in modern C++ (C ++14/17) using [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. I wrote it while studying book
[Practical mathematical optimization: An Introduction to Basic Optimization Theory and Classical and New Gradient-Based Algorithms](https://www.springer.com/gp/book/9780387298245)
by Jan A.Snyman.
## Content
### Unconstrained optimization
* One-Dimensional Line Search Methods (Golden Section, Fibonacci, Brent, Powell)
* Steepest Descent Method
* Conjugate Gradient Methods (Fletcher-Powell, Polak-Ribiere, Hestenes-Stiefel)
* Quasi-Newton Methods (Davidon-Fletcher-Powell, Broyden-Fletcher-Goldfarb-Shanno)
* Zero-Order Methods (Nelder-Mead, Powell's Conjugate Method)
## Usage
See uni_tests folder