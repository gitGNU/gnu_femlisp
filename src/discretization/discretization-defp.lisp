;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; discretization-defp.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are
;;; met:
;;; 
;;; 1. Redistributions of source code must retain the above copyright
;;; notice, this list of conditions and the following disclaimer.
;;; 
;;; 2. Redistributions in binary form must reproduce the above copyright
;;; notice, this list of conditions and the following disclaimer in the
;;; documentation and/or other materials provided with the distribution.
;;; 
;;; THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
;;; WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defpackage "FL.DISCRETIZATION"
  (:nicknames "DISC")
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.MESH" "FL.PROBLEM" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.DEBUG" "FL.ITERATION")
  
  (:export				; discretization.lisp
   "<DISCRETIZATION>" "DISCRETIZATION-ORDER" "DISCRETIZE"
   "SELECT-DISCRETIZATION")
  
  (:export				; fe.lisp
   "<DOF>" "DOF-INDEX" "DOF-SUBCELL-INDEX" "DOF-IN-VBLOCK-INDEX"
   "DOF-VBLOCK-LENGTH"
   "DOF-FUNCTIONAL" "DOF-COORD" "DOF-GCOORD" "INTERIOR-DOF?" "DOF-COMPONENT"
   "ASSEMBLE-INTERIOR" "ESSENTIAL-BOUNDARY-CONSTRAINTS" "CONSTRAINTS-ON-LEVEL"
   "ELIMINATE-CONSTRAINTS"
   "<FE>" "FE-DOFS" "FE-BASIS" "INTEGRATION-POINTS"
   "MAKE-LOCAL-VEC" "MAKE-LOCAL-MAT" "FE-CELL-GEOMETRY"
   "NR-OF-DOFS" "NR-OF-INNER-DOFS"
   "SUBCELL-NDOFS" "SUBCELL-INDICES" "INNER-DOF-INDICES"
   "INTERPOLATE-ON-REFCELL"
   "<VECTOR-FE>" "COMPONENTS" "NR-OF-COMPONENTS"
   "SUBCELL-OFFSETS" "LOCAL-OFFSET"
   "<FE-DISCRETIZATION>" "<STANDARD-FE-DISCRETIZATION>"
   "<SCALAR-FE-DISCRETIZATION>" "<VECTOR-FE-DISCRETIZATION>"
   "GET-FE" "CELL->FE" "QUADRATURE-RULE" "IP-VALUES" "IP-GRADIENTS"
   "LOCAL-IMATRIX" "LOCAL-PMATRIX" "LOCAL-TRANSFER-MATRIX"
   "COMPONENT")
  
  (:export				; fedisc.lisp
   "ASSEMBLE-CONSTRAINTS" "FE-DISCRETIZE" "DISCRETIZE-LOCALLY" "DISCRETIZE-GLOBALLY")

  (:export				; feeval.lisp
   "FE-LOCAL-VALUE" "FE-VALUE"
   "FE-LOCAL-GRADIENT" "FE-GRADIENT"
   "FE-INTEGRATE" "CELL-INTEGRATE"
   "FE-EXTREME-VALUES")

  (:export				; lagrange.lisp
   "LAGRANGE-FE" "LAGRANGE-VECTOR-FE" "LAGRANGE-MAPPING")

  (:export				; quadrature.lisp
   "<IP>" "IP-WEIGHT" "IP-COORDS" "GAUSS-RULE"
   "GAUSS-LOBATTO-POINTS-ON-UNIT-INTERVAL")
   
  (:export				; sparseif
   "GET-LOCAL-FROM-GLOBAL-VEC"
   "SET-GLOBAL-TO-LOCAL-VEC" "INCREMENT-GLOBAL-BY-LOCAL-VEC"
   "CELL->MATRIX-VALUE-BLOCKS" "GET-LOCAL-FROM-GLOBAL-MAT"
   "SET-GLOBAL-TO-LOCAL-MAT" "INCREMENT-GLOBAL-BY-LOCAL-MAT")

   (:export				; sparseas.lisp
    "<ANSATZ-SPACE>" "MAKE-FE-ANSATZ-SPACE" "PROBLEM"
    "CELL-KEY" "REPRESENTATIVE"
    "<ANSATZ-SPACE-VECTOR>" "MAKE-ANSATZ-SPACE-VECTOR"
    "RANDOM-ANSATZ-SPACE-VECTOR"
    "<ANSATZ-SPACE-MORPHISM>" "MAKE-ANSATZ-SPACE-MORPHISM"
    "<ANSATZ-SPACE-AUTOMORPHISM>" "MAKE-ANSATZ-SPACE-AUTOMORPHISM"
    "ANSATZ-SPACE" "FE-CLASS" "MESH" "HIERARCHICAL-MESH"
    "SURFACE-CELLS" "EXTRACT-LEVEL"
    "INTERPOLATION-MATRIX" "CONSTRAINED-INTERPOLATION-MATRIX" "PROJECTION-MATRIX"
    "TRANSFER-MATRIX" "SORT-KEYS" "DECOMPOSE")
   
   (:documentation "The @package{FL.DISCRETIZATION} package defines
@class{<discretization>} as an abstract class and
@class{<fe-discretization>} as a concrete derived class.

The key for local assembly is given by the generic function
@function{get-fe}, which yields a suitable finite element for a given cell.
The value of @function{get-fe} is a class @class{<fe>} for scalar problems
or @class{<vector-fe>} for vector-valued problems which contains
information on base functions and node functionals.  Another generic
function @function{quadrature-rule} computes memoized quadrature rules for
those finite elements.

Obviously, also non-standard finite element discretizations like hp-methods
would fit into this scheme.  The key for local assembly is given by the
generic function @function{get-fe}, which yields a suitable finite element
for a given cell.  The value of @function{get-fe} is a class @class{<fe>}
for scalar problems or @class{<vector-fe>} for vector-valued problems which
contains information on base functions and node functionals.  Another
generic function @function{quadrature-rule} computes memoized quadrature
rules for those finite elements.

The file @path{lagrange.lisp} provides Lagrange finite elements of
arbitrary order.  The evaluation points for the node functionals may be
chosen either uniformly distributed or at the Gauss-Lobatto points.  The
latter choice of points yields better stability properties but is
restricted to cube meshes.  Also functions for constructing cell mappings
by pointwise evaluation of the domain boundary are provided here, which may
be used to construct isoparametric domain approximations.

In the file @path{fedisc.lisp}, the function @function{fe-discretize} is
defined.  This function performs the standard steps for finite element
discretization: interior assembly, assembly and elimination of hanging-node
and essential-boundary constraints.  It works on a blackboard as explained
in Section @ref{Blackboards} and can reuse already available matrix-vector
structure.  There is a somewhat less flexible interface provided by the
funtion @function{discretize-globally} which calls
@function{fe-discretize}.

In the files @path{cdr-fe.lisp}, @path{elasticity-fe.lisp} and
@path{navier-stokes.lisp} one can find methods for local assembly for the
different problems.  They are defined in own packages which use both the
package @package{FL.DISCRETIZATION} and the package for the particular
problem."))


