(in-package :ddo-femlisp)

#+(or)
(defun convert-to-full-tensor (arr/mat)
  (lret* ((dim1(array-dimension arr/mat 0))
             (tensi (make-real-tensor (make-array 4 :initial-element dim))))
       (for-each-entry-and-key
        (lambda (x &rest indices)
       (for-each-entry-and-key
        (lambda (value &rest indices2)
          (apply '(setf tensor-ref) value tensi (append indices indices2)))
        x))
     arr/mat)))

#+(or)
(defun copy-full-tensor-to-arr/mat (tensor result)
  (for-each-entry-and-key
   (lambda (mat i j)
     (for-each-key
      (lambda (k l)
        (setf (mref mat k l) (tensor-ref tensor i j k l)))
      mat))
   result)
  result)

(add-hook 'fl.application::effective-tensor
          'distribute-if-necessary
          (lambda (result as)
            (when (get-property as :distributed-p)
              (loop for object in (etypecase result
                                    (standard-matrix (list result))
                                    (array (entries result)))
                    do
                       (make-distributed-object
                        object (all-processors) '(store)))
              (let ((*synchronization-merger* #'addition-merger))
                (synchronize))
            result)))


