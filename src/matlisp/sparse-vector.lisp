
(defclass container ()
  ((entries)
   (free-indices)
   (indexer :type indexer)))

(defclass indexer ()
  (key->index ...))

(defclass typed-container ()
  ((entries)
   (free-indices)
   (indexer)))


(loop with i = 0 and free-indices = free-indices
      and free-index = (car free-indices)
      
1. Array
2. Hash-Tabelle 