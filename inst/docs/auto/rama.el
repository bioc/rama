(TeX-add-style-hook "rama"
 (function
  (lambda ()
    (LaTeX-add-bibliographies
     "robust_estimation")
    (LaTeX-add-labels
     "tab:data_format")
    (TeX-run-style-hooks
     "Sweave"
     "chicago"
     "hyperref"
     "graphicx"
     "amsmath"
     "epsfig"
     "psfig"
     "fullpage"
     "latex2e"
     "art11"
     "article"
     "11pt"))))

