(in-package cl-user)

(defpackage :cl-numerics-utils
  (:use #:cl #:cl-linear-algebra)
  (:export #:num= #:*num=-tolerance*
           #:ewt #:ewt-rtol #:ewt-atol #:make-ewt #:ewt-p
           #:ewt-vector
           #:rms-norm #:rms-norm-diff
           #:horner-rule))
