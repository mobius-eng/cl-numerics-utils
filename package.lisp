(in-package cl-user)

(defpackage :cl-numerics-utils
  (:use #:cl #:cl-linear-algebra)
  (:export #:plus-infinite-p #:minus-infinite-p
           #:+double-float-plus-infinity+
           #:+double-float-minus-infinity+
           #:+single-float-plus-infinity+
           #:+single-float-negative-infinity+
           #:num= #:*num=-tolerance*
           #:ewt #:ewt-rtol #:ewt-atol #:make-ewt #:ewt-p
           #:ewt-vector
           #:rms-norm #:rms-norm-diff
           #:horner-rule))
