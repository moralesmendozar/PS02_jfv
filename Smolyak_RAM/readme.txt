Code for Computing Dynamic Heterogeneous-Agent Economies: Tracking the Distribution

Available at https://sites.google.com/site/greygordon

This Matlab code performs Smolyak Interpolation in three steps.  
(1) smolyakapprox_step1 : this provides collocation points and a structure which is basically a template for constructing the polynomial.
(2) smolyakapprox_step2 : given the function values and the structure from (1), this constructs the polynomial coefficients placing them in a new structure.
(3) smolyakapprox_step3 : given the structure from (2), this evaluates the polynomial at arbitrary points in the state space.

Note: I plan on coming back and providing better documentation as well as the code to reproduce these results as soon as the job market is over.  For now I hope this will do.

All this code is free to modify under the GNU license below

    Code for Computing Dynamic Heterogeneous-Agent Economies: Tracking the Distribution
    Copyright (C) 2011  Grey Gordon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
