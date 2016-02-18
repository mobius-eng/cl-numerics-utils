(asdf:defsystem #:cl-numerics-utils
  :author "Alexey Cherkaev (mobius-eng)"
  :license "LGPLv3"
  :description "Various utils for numerical methods"
  :serial t
  :depends-on (:cl-linear-algebra)
  :components ((:file "package")
               (:file "cl-numerics-utils")))
