# FFD.jl
A Julia repository containing code for Free Form Deformation using B-splines.
The repository has the following structure

* **FFD.jl**
  * **src**
    * bounding_box.jl
    * b-splines.jl
    * control_point.jl
    * evaluations.jl
    * extra.jl
    * knot.jl
    * mapping_functions.jl
    * mapping.jl
    * param_functions.jl
    * plot.jl
    * span.jl
    * startup.jl
  * **test**
    * runtest.jl
    * test_b-splines.jl
    * test_linearFFD.jl
  * README.md

The code was developed from the following references

1.  Les Piegl and Wayne Tiller. 1997. *The NURBS Book (2nd Ed.)*.
Springer-Verlag New York, Inc., New York, NY, USA.
2.  Carl de Boor. 1978. *A Prcatical Guide to Splines*. Springer-Verlag
New York, Inc., New York, NY, USA.
3.  Thomas W. Sederberg and Scott R. Parry. 1986. Free-form deformation of solid
geometric models. In *Proceedings of the 13th annual conference on Computer
graphics and interactive techniques* (SIGGRAPH '86), David C. Evans and Russell
J. Athay (Eds.). ACM, New York, NY, USA, 151-160.
DOI=http://dx.doi.org/10.1145/15922.15903
4. Prof. Hicken's Mapping_Mod.f90
5. University of Michigan Free-form Deformation Toolbox
