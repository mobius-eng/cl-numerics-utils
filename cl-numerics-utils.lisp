(in-package cl-numerics-utils)
;; * Common utils for numerical methods
;; ** Infinities

(defgeneric plus-infinite-p (number))
(defgeneric minus-infinite-p (number))

(defmethod plus-infinite-p ((number number)) nil)
(defmethod minus-infinite-p ((number number)) nil)


(defmethod plus-infinite-p ((number double-float))
  (> number most-positive-double-float))

(defmethod minus-infinite-p ((number double-float))
  (< number most-negative-double-float))


(defmethod plus-infinite-p ((number single-float))
  (> number most-positive-single-float))

(defmethod minus-infinite-p ((number single-float))
  (< number most-negative-single-float))


(defconstant +double-float-plus-infinity+
  #+sbcl sb-ext:double-float-positive-infinity
  #+ccl 1d++0
  #+ecl ext:double-float-positive-infinity)

(defconstant +double-float-minus-infinity+
  #+sbcl sb-ext:double-float-negative-infinity
  #+ccl -1d++0
  #+ecl ext:double-float-negative-infinity)

(defconstant +single-float-plus-infinity+
  #+sbcl sb-ext:single-float-positive-infinity
  #+ccl 1e++0
  #+ecl ext:single-float-positive-infinity)

(defconstant +single-float-negative-infinity+
  #+sbcl sb-ext:single-float-negative-infinity
  #+ccl -1e++0
  #+ecl ext:single-float-negative-infinity)


;; ** Generalized equality: borrowed from CL-NUM-UTILS
(defgeneric num= (x y &optional tolerance)
  (:documentation "Test if X equals Y within the precision TOLERANCE.
    How TOLERANCE is defined or specified may vary depending on types of X and Y"))

(defvar *num=-tolerance* (sqrt double-float-epsilon)
  "Default tolerance for double-float numbers")

(defmethod num= ((x number) (y number) &optional (tolerance *num=-tolerance*))
  (< (abs (- x y)) tolerance))

(defclass ewt ()
  ((rtol :initarg :rtol
         :reader ewt-rtol
         :documentation "Relative tolerance")
   (atol :initarg :atol
         :reader ewt-atol
         :documentation "Absolute tolerance")
   (result :initarg :result
           :reader ewt-result
           :documentation "Internal storage for result"))
  (:documentation "Error weight vector"))

(defmethod initialize-instance :after ((obj ewt) &key)
  (with-slots (rtol result) obj
    (setf result (make-vector (length rtol) 'double-float 0d0))))

(defgeneric ewt (rtol atol &key &allow-other-keys)
  (:documentation "Various ways of constructing EWT"))

(defmethod ewt ((rtol vector) (atol vector) &key)
  (check-vector-lengths rtol atol)
  (make-instance 'ewt :atol atol :rtol rtol))

(defmethod ewt ((rtol vector) (atol number) &key)
  (make-instance 'ewt
    :atol (make-vector (length rtol) (type-of atol) atol)
    :rtol rtol))

(defmethod ewt ((rtol number) (atol vector) &key)
  (make-instance 'ewt
    :atol atol
    :rtol (make-vector (length atol) (type-of rtol) rtol)))

(defmethod ewt ((rtol vector) atol &key)
  (ewt rtol *num=-tolerance*))

(defmethod ewt (rtol (atol vector) &key)
  (ewt *num=-tolerance* atol))

(defmethod ewt ((rtol number) (atol number) &key size)
  (ewt (make-vector size 'double-float (coerce rtol 'double-float))
       (make-vector size 'double-float (coerce atol 'double-float))))

(defmethod ewt (rtol atol &key size)
  (ewt *num=-tolerance* *num=-tolerance* :size size))

(defun ewt-p (obj) (eq (type-of obj) 'ewt))

(defun ewt-vector (ewt y)
  "Constructs error weight vector from EWT tolerances and vector Y."
  (let ((z (ewt-result ewt)))
    (vector-elementwise z #'(lambda (p q) (abs (* p q))) (ewt-rtol ewt) y)
    (add-vectors z :vectors (list (ewt-atol ewt)) :multipliers (list 1))
    z))

(defun rms-norm (vector ewt-vector)
  "Weighted root-mean-square norm

          -----------------------
        /    ---  |    v    | 2
       /  1  \\    |     i   |
      /   -  /    | ------- |
  \\  /    N  ---  |  EWT    |
   \\/         N   |     i   |

"
  (let ((N (length vector))
        (result 0))
    (iter-vector (((x vector) (y ewt-vector)) (i) (sqrt (* (/ N) result)))
      (incf result (expt (/ (x i) (y i)) 2)))))

(defun rms-norm-diff (vector1 vector2 ewt-vector)
  "Weighted root-mean-square norm of difference of two vectors"
  (let ((N (length vector1))
        (result 0))
    (iter-vector (((x1 vector1) (x2 vector2) (y ewt-vector)) (i) (sqrt (* (/ N) result)))
      (incf result (expt (/ (- (x1 i) (x2 i)) (y i)) 2)))))

(defmethod num= ((x vector) (y vector) &optional tolerance)
  (let ((tolerance (cond ((or (null tolerance) (numberp tolerance) (vectorp tolerance))
                          (ewt tolerance tolerance :size (length x)))
                         ((ewt-p tolerance) tolerance)
                         ((consp tolerance)
                          (let ((rtol (car tolerance))
                                (atol (if (consp (cdr tolerance))
                                          (cadr tolerance)
                                          (cdr tolerance))))
                            (ewt rtol atol :size (length x))))
                         (t (error "Unknown tolerance type: ~A" tolerance)))))
    (<= (rms-norm-diff x y (ewt-vector tolerance x)) 1d0)))


(defgeneric horner-rule (polynomial number)
  (:documentation "Horner's rule of polynomial evaluation"))

(defmethod horner-rule ((polynomial cons) number)
  "Polynomial coefficients must be stored from heighest to lowerest powers"
  (let ((result (first polynomial)))
    (dolist (coeff (rest polynomial) result)
      (setf result (+ coeff (* number result))))))

(defmethod horner-rule ((polynomial vector) number)
  "Polynomial coefficients must be stored from lowest to heighest powers"
  (let* ((length (length polynomial))
         (result (aref polynomial (1- length))))
    (loop for i from (- length 2) downto 0
       do (setf result (+ (aref polynomial i) (* number result)))
       finally (return result))))

