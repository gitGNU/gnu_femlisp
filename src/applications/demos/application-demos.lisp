;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; application-demos.lisp
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

(in-package :fl.application)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refinement demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *refinement-demos*
  (make-demo
   :name "Refine"
   :short "Shows regular refinements of geometric bodies."
   :long
   "The refinement is done using Freudenthal's algorithm for
simplex refinement in a form generalized to tensor product
cells."))

(adjoin-demo *refinement-demos* *demo-root*)

(defvar *discretization-demo*
  (make-demo
   :name "Discretization"
   :short "Demonstrate discretization aspects"))
(adjoin-demo *discretization-demo* *demo-root*)

(defvar *solver-demo*
  (make-demo
   :name "Solvers"
   :short "Demonstrate solver aspects"))
(adjoin-demo *solver-demo* *demo-root*)

(defvar *equation-demo*
  (make-demo
   :name "Equations"
   :short "Solve several types of pdes"))
(adjoin-demo *equation-demo* *demo-root*)

(defvar *adaptivity-demo*
  (make-demo
   :name "Adaptivity"
   :short "Adaptive solution"))
(adjoin-demo *adaptivity-demo* *demo-root*)

(defvar *effective-coeffs-demo*
  (make-demo
   :name "Homogenization"
   :short "Computing effective coefficients"
   :long
   "In homogenization for periodic media the computation of
effective constants is done via solution of a cell problem in a
representative cell."))
(adjoin-demo *effective-coeffs-demo* *demo-root*)

(defvar *interior-coeffs-demo*
  (make-demo
   :name "interior-coeffs"
   :short "Computing interior coefficients"
   :long
   "The effect of periodical oscillations in coefficients of a
partial differential equation can be replaced up to some error
by an averaged equation.  Computing the coefficients appearing
in these effective equations usually involves solving a cell
problem on a representative cell and averaging the obtained
solution."))
(adjoin-demo *interior-coeffs-demo* *effective-coeffs-demo*)

(defvar *boundary-coeffs-demo*
  (make-demo
   :name "boundary-coeffs"
   :short "Computing boundary law coefficients"
   :long
   "The effect of a periodically oscillating boundary can be
replaced up to some error by including correction terms into the
boundary law.  Computing the coefficients appearing in these
usually involves solving a problem on a semi-infinite
periodicity cell and averaging the obtained solution."))
(adjoin-demo *boundary-coeffs-demo* *effective-coeffs-demo*)

(defvar *articles-demo*
  (make-demo
   :name "Articles"
   :short "Calculations connected to articles"
   :long
   "The main goal of Femlisp is education and research.  In this
part of the demonstrations, we have collected computations which
were done in Femlisp related to some research articles."))
(adjoin-demo *articles-demo* *demo-root*)

(defvar *books-demo*
  (make-demo
   :name "Books"
   :short "Collection of demos for books"
   :long
   "The main goal of Femlisp is education and research.  In this
part of the demonstrations, we have collected computations which
are related to some books."))
(adjoin-demo *books-demo* *demo-root*)

(defvar *courses-demo*
  (make-demo
   :name "Courses"
   :short "Collection of demos for courses"
   :long
   "A collection of demos for courses about different
subjects."))
(adjoin-demo *courses-demo* *demo-root*)

(defvar *talks-demo*
  (make-demo
   :name "Talks"
   :short "Collection of demos for talks"
   :long
   "A collection of demos for talks about different subjects."))
(adjoin-demo *talks-demo* *demo-root*)




