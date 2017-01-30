(in-package :femlisp-ddo)

(defun diagonal-asa (asa)
  (lret ((result (make-analog asa)))
    (for-each-row-key
     (lambda (rk)
       (copy! (mref asa rk rk) (mref result rk rk)))
     asa)
    (aso-make-distributed result)))

(defun make-consistent (vec)
  (when (get-property vec :distributed-p)
    (dovec (entry vec)
      (when (distributed-p entry)
        (insert-into-changed entry)))
    (let ((*synchronization-merger* #'addition-merger))
      (synchronize 'make-consistent vec)))
  vec)
  
(defun jacobi-correction (sol omega corr diag res)
  (copy! res corr)
  (dorows (rk diag)
    (gesv! (mref diag rk rk) (vref corr rk)))
  (make-consistent corr)
  (axpy! omega corr sol))

(defclass <distributed-jacobi> (fl.iteration::<linear-iteration>)
  ())

(defmethod fl.iteration::make-iterator
    ((jac <distributed-jacobi>) (mat <ansatz-space-automorphism>))
  (let (corr
        (diag (diagonal-asa mat)))
    (make-instance
     'fl.iteration::<iterator>
     :matrix mat
     :residual-before t
     :initialize nil
     :iterate
     #'(lambda (x b r)
         (declare (ignore b))
         (unless corr
           (setf corr (make-domain-vector-for mat)))
         (jacobi-correction x (slot-value jac 'fl.iteration::damp) corr diag r))
     :residual-after nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <distributed-cg> (fl.iteration::<linear-iteration>)
  ((preconditioner :initform nil :initarg :preconditioner)
   (restart-cycle :reader restart-cycle :initform nil :initarg :restart-cycle))
  (:documentation "Preconditioned conjugate gradient iteration"))

(defmethod fl.iteration::make-iterator ((cg <distributed-cg>) mat)
  "Standard method for the preconditioned conjugate-gradient iteration."
  (declare (optimize debug))
  (let* ((precond (aand (slot-value cg 'preconditioner)
			(fl.iteration::make-iterator it mat))))
    (when (and precond (slot-value precond 'fl.iteration::residual-after))
      (error "This preconditioner does not work, because the application
here wants to keep residual and rhs intact."))
    (let ((p (make-domain-vector-for mat))
	  (a (make-image-vector-for mat))
	  (w (make-image-vector-for mat))
	  (q (make-domain-vector-for mat))
	  (alpha 0.0)
          (count 0))
      (assert (and p a w))
      (with-slots (fl.iteration::initialize fl.iteration::iterate) precond
        (flet ((restart (r)
                 (cond
                   (precond
                    (copy! r w)
                    (when fl.iteration::initialize
                      (funcall fl.iteration::initialize p w w))
                    (funcall fl.iteration::iterate p w w))
                   (t (copy! r p)
                      (make-consistent p)))
                 (setq alpha (dot p r))))
          (make-instance
           'fl.iteration::<iterator>
           :matrix mat
           :residual-before t
           :initialize
           #'(lambda (x b r)
               (declare (ignore x b))
               (restart r))
           :iterate
           #'(lambda (x b r)
               (declare (ignore b))
               (when (aand (plusp count)
                           (restart-cycle cg)
                           (zerop (mod count it)))
                 (restart r))
               (unless (zerop alpha)
                 (gemm! 1.0 mat p 0.0 a)
                 (let* ((beta (dot a p))
                        (lam (/ alpha beta)))
                   (axpy! lam p x)
                   (axpy! (- lam) a r)
                   (let ((q (cond (precond
                                   (copy! r w) (x<-0 q)
                                   (funcall fl.iteration::iterate q w w)
                                   q)
                                  (t (copy! r q)
                                     (make-consistent q)))))
                     (let ((new-alpha (dot q r)))
                       (scal! (/ new-alpha alpha) p)
                       (m+! q p)
                       (setq alpha new-alpha))
                     )))
               ;; restart procedure
               (unless (aand (restart-cycle cg)
                             (zerop (mod (incf count) it)))
                 (list :residual-after t)))
           :residual-after t))))))

(defun ensure-consistent-copy (x)
  (unless (get-property x :distributed-p)
    (when (get-property (ansatz-space x) :distributed-p)
      (or (get-property x 'consistent-copy)
          (lret ((xc (copy x)))
            (aso-make-distributed xc)
            (setf (get-property x 'consistent-copy) xc))))))

(defmethod l2-norm ((x <ansatz-space-vector>))
  (assert (not (get-property x :distributed-p)))
  (if (get-property (ansatz-space x) :distributed-p)
      (let ((xc (or (ensure-consistent-copy x) x)))
        (copy! x xc)
        (make-consistent xc)
        (sqrt (dot x xc)))
      (call-next-method)))

(defmethod dot :around ((x <ansatz-space-vector>) (y <ansatz-space-vector>))
  (let ((result (call-next-method)))
    (cond ((or (get-property x :distributed-p)
               (get-property y :distributed-p))
           (let ((comm (make-distributed-object
                        (ensure-matlisp result) (all-processors) '(store))))
             (when (and (get-property x :distributed-p)
                        (get-property y :distributed-p))
               ;; unfortunately, we need a correction step in this case.
               ;; ideally, an algorithm should be written such that this kind
               ;; of dot product is not necessary
               (for-each-entry-and-key
                (lambda (xc i)
                  (when (distributed-p xc)
                    (unless (masterp xc)
                      (decf (vref comm 0) (dot xc (vref y i))))))
                x))
             (let ((*synchronization-merger* #'addition-merger))
               (synchronize 'dot))
             (vref comm 0)))
          (t result))))
