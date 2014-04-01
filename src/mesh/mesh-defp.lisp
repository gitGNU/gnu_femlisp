;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; mesh-defp.lisp
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

(defpackage "FL.MESH"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION")
  (:export
   ;; cell.lisp
   "<CELL>" "CELLP" "BOUNDARY"
   "<MAPPED-CELL>" "MAPPED-P" "MAPPING" "MAPPED-CELL-CLASS"
   "DIMENSION" "EMBEDDED-DIMENSION"
   "NR-OF-SIDES" "NR-OF-VERTICES" "NR-OF-SUBCELLS"
   "FACTOR-SIMPLICES" "FACTOR-DIMENSIONS" "REFERENCE-CELL" "REFERENCE-CELL-P"
   "SUBCELLS" "VERTICES" "CORNERS" "DIAMETER"
   "CELL-MAPPING"
   "L2G" "L2DG"
   "LOCAL->GLOBAL" "LOCAL->DGLOBAL"
   "MULTIPLE-LOCAL->GLOBAL" "MULTIPLE-LOCAL->DGLOBAL"
   "G2L" "GLOBAL->LOCAL" "GLOBAL->EMBEDDED-LOCAL"
   "INSIDE-CELL?" "LOCAL-COORDINATES-OF-MIDPOINT" "MIDPOINT" "ORIGIN"
   "<VERTEX>" "MAKE-VERTEX" "VERTEX-POSITION" "*REFERENCE-VERTEX*"
   "VERTEX?" "VERTEX-P" "PRINT-VERTEX"
   "MAKE-CELL-FROM-VERTICES" "MAKE-CELL-FROM-CORNERS"

   ;; simplex.lisp
   "<SIMPLEX>" "SIMPLEX-CLASS" "SIMPLEX-P" "MAKE-SIMPLEX" "MAKE-LINE"
   "*UNIT-INTERVAL*" "*UNIT-TRIANGLE*" "*UNIT-TETRAHEDRON*"
   "N-SIMPLEX"
   ;; product-cell.lisp
   "<PRODUCT-CELL>" "*UNIT-QUADRANGLE*" "*UNIT-CUBE*" "N-CUBE" "CUBE-P"
   "CELL->CUBE" "ENSURE-SIMPLEX-PRODUCT" "CARTESIAN-PRODUCT"
   
   ;; skeleton.lisp
   "<SKELETON>" "SKELETON" "ETABLES" "ETABLE" "SKEL-EMPTY-P"
   "SKEL-REF" "GET-CELL-PROPERTY" "MARK-SKELETON"
   "INSERT-CELL!" "INSERT-CELLS!"
   "NR-OF-CELLS" 
   "SKEL-FOR-EACH" "DOSKEL" "SKEL-MAP" "SKEL-FULL-MAP"
   "FIND-CELLS" "FIND-CELL" "FIND-CELL-FROM-POSITION" "FIND-CELL-FROM-CORNERS"
   "CELLS-OF-DIM" "ETABLE-OF-HIGHEST-DIM"
   "CELLS-OF-HIGHEST-DIM" "SURFACE-CELLS-OF-HIGHEST-DIM"
   "FOR-EACH-CELL-OF-HIGHEST-DIMENSION-ON-SURFACE"
   "CHECK"
   "MEMBER-OF-SKELETON?" "SKELETON-BOUNDARY"
   "SUBSKELETON"

   ;; identify
   "CELL-IDENTIFICATION" "REPRESENTATIVE"
   "IDENTIFIED-P" "IDENTIFIED-CELLS" "KEY-IS-SUBCELL-P"
   "BOUNDARY-IDENTIFICATIONS" "COMBINE-IDENTIFICATIONS"
   "ITERATE-IDENTIFICATIONS" "IDENTIFY" 
   "IDENTIFY-UNIT-CELL-FACES"
   
   ;; skeleton-build.lisp
   "SKELETON-DISJOINT-UNION" "SKELETON-WITHOUT-CELL" "SKEL-ADD!" "COPY-SKELETON"
   "TRANSFORMED-SKELETON" "LINEARLY-TRANSFORMED-SKELETON" "SHIFT-SKELETON"
   "TELESCOPE" "REFINED-SKELETON-COPY" "STRUCTURED-SKELETON"
   
   ;; refine.lisp
   "REFINEMENT" "REFINEMENT-RULE" "GET-REFINEMENT-RULE"
   "CHILDREN" "PARENT" "REFINED-P"
   "REFINE-INFO" "SKELETON-REFINEMENT" "INITIALIZE-REFINED-SKELETON"
   "INVERT-REFINEMENT"
   "SUBCELL-CHILDREN" "REFCELL-CHILDREN" "INNER-REFCELL-CHILDREN"
   "REFCELL-REFINEMENT-SKELETON"

   ;; domain.lisp
   "<DOMAIN>" "DOMAIN-BOUNDARY"
   "MAKE-CLASSIFIER" "DOMAIN-CHARACTERISTICS" "EXTENSIBLE-P"
   "MAKE-DOMAIN"
   "PATCH" "PATCH-OF-CELL" "PATCH-CLASSIFICATION" "TEST-CONDITION" "FIND-PATCH"
   "DIMENSION-OF-PART"
   "N-SIMPLEX-DOMAIN" "SIMPLEX-PRODUCT-DOMAIN" "N-CUBE-DOMAIN" "ENSURE-DOMAIN"
   "BOX-DOMAIN" "N-CELL-DOMAIN" "N-BALL-DOMAIN"
   "TRIANGLE-DOMAIN" "L-DOMAIN"

   ;; mesh.lisp
   "<MESH>" "DOMAIN"
   "REFINE" "<HIERARCHICAL-MESH>" "MESHES" "REFINEMENTS" "ACCUMULATED-MESH"
   "NR-OF-LEVELS" "CELLS-ON-LEVEL" "TOP-LEVEL" "BOTTOM-LEVEL-CELLS" "TOP-LEVEL-CELLS"
   "SURFACE-CELLS-OF-DIM" "NR-OF-SURFACE-CELLS"
   "LEVEL-OF-CELL" "MESHSIZE"
   "MESH->HIERARCHICAL-MESH" "REFINE-HIERARCHICAL-MESH"
   "REFINEMENT-INTERFACE" "HIERARCHICALLY-ORDERED-CELLS"
   
   ;; meshgen.lisp
   "MAKE-MESH-FROM" "COPY-MESH" "UNIFORMLY-REFINED-MESH"
   "INSERT-CELL-FROM-CORNERS" "SPECIAL-MESH-ON-BOX-DOMAIN"
   "UNIFORM-MESH-ON-BOX-DOMAIN"
   
   "MAKE-HIERARCHICAL-MESH-FROM-DOMAIN" "MAKE-HIERARCHICAL-MESH-FROM"
   "UNIFORMLY-REFINED-HIERARCHICAL-MESH"
   "COMPARE-LEXICOGRAPHICALLY" "SORT-LEXICOGRAPHICALLY"
   "TRIANGULIZE"

   ;; triangulate.lisp
   "<BOUNDARY-CELL>" "TRIANGULATE"
   
   ;; extend.lisp
   "EXTENSION" "EXTEND" "STANDARD-EXTENDER" "CUBE-EXTENDER"
   )
  (:documentation "This module contains the definitions of meshes and
routines for mesh management.  The meshes allowed in @femlisp{} are more
general than those of most other software for solving PDEs.  In @femlisp{},
both mesh, domain and problem definitions are defined over an underlying
abstraction, the so-called @class{<skeleton>}.  A @class{<skeleton>}
captures the mathematical idea of a \"cell complex\" which builds a
topological space by mapping from standard cells @class{<cell>}.  Now, a
@class{<skeleton>} can be seen as mapping the cells of such a cell complex
to arbitrary values.  Then, a @class{<domain>} is a @class{<skeleton>}
where each cell (which we call \"patch\" in this case) is mapped to
geometric properties, and a @class{<mesh>} is a @class{<skeleton>} where
each cell is mapped to the patch to which it belongs.

The basic entities are the class @class{<cell>}, the subclass
@class{<simplex>} which in turn contains subclasses for arbitrarily
dimensional simplices generated on demand, and the subclass
@class{<product-cell>} containing arbitrary products of simplices, e.g. square
or cube.

Meshes can be refined either uniformly or locally using the Freudenthal
algorithm as presented in @cite{JBey_2000a} and generalized to product
elements.  When local refinement is used, hanging nodes may occur.  In
contrast to most other finite element software, in @femlisp{} the
difference of refinement levels of adjacent cells may be arbitrarily large.
Up to now, anisotropic refinement of product cells has not yet been
implemented."))
