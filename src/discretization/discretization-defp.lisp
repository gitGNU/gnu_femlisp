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
	"FL.DEBUG")
  
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
    "<ANSATZ-SPACE-MORPHISM>" "MAKE-ANSATZ-SPACE-MORPHISM"
    "<ANSATZ-SPACE-AUTOMORPHISM>" "MAKE-ANSATZ-SPACE-AUTOMORPHISM"
    "ANSATZ-SPACE" "FE-CLASS" "MESH" "HIERARCHICAL-MESH"
    "SURFACE-CELLS" "EXTRACT-LEVEL"
    "INTERPOLATION-MATRIX" "CONSTRAINED-INTERPOLATION-MATRIX" "PROJECTION-MATRIX"
    "TRANSFER-MATRIX" "SORT-KEYS" "DECOMPOSE")
   
   (:documentation "The @code{FL.DISCRETIZATION} package defines
@code{<discretization>} as an abstract class and @code{<fe-discretization>}
as a concrete derived class.  A generic function @code{get-fe} is used for
associating a cell with a finite element @code{<fe>}, which is a data
structure containing information about base functions and dual functionals
on the corresponding cell.  Lagrange finite elements of arbitrary order are
implemented as a special instance of @code{<fe-discretization>}.  Note that
other discretizations as finite differences or finite volumes could easily
be incorporated as well."))

