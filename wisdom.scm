(define h-plausible 0.000001)

(define (linear-stability map x0 y0 eps result)
  (map x0 y0 
       (lambda (nx0 ny0)
	 (if (and (< (square (- nx0 x0)) (square eps))
		  (< (square (- ny0 y0)) (square eps)))
	     (let ((a11 ((richardson-derivative (lambda (x) (map x y0 (lambda (nx ny) nx) (lambda () (result 'error)))) eps h-plausible) x0))
		   (a22 ((richardson-derivative (lambda (y) (map x0 y (lambda (nx ny) ny) (lambda () (result 'error)))) eps h-plausible) y0))
		   (a12 ((richardson-derivative (lambda (x) (map x y0 (lambda (nx ny) ny) (lambda () (result 'error)))) eps h-plausible) x0))
		   (a21 ((richardson-derivative (lambda (y) (map x0 y (lambda (nx ny) nx) (lambda () (result 'error)))) eps h-plausible) y0)))
	       (if (> (- (square (+ a11 a22)) 4) 0) (result 'unstable) (result 'stable)))
	     (error "not a fixed point" x0 y0 nx0 ny0)))
       (lambda () (result 'error))))

(define (read-nth nth)
  (if (= nth 0)
      (read)
      (begin (read)
	     (read-nth (- nth 1)))))

(define (HEMdoit6 name a1/Re inc Ie H0+H1/L0 eps cont)
  (HEMconstants name a1/Re inc Ie H0+H1/L0)
  (let ((map (HEMmap name 0.1)))
    (let ((Ls/L0 (with-input-from-file (string-append name ".constants") (lambda () (read-nth 8)))))
      (linear-stability map 0. 0. eps (cont Ls/L0)))))

(define (HEMdoit-range6 a1/Re H0+H1/L0 eps imin imax di Iemin Iemax dIe)
  (let ((name (string-append "tmp" (number->string (random 100000)))))
    (let iloop ((inc imin))
      (if (< inc imax)
	  (let Ieloop ((Ie Iemin))
	    (if (< Ie Iemax)
		(begin
		  (HEMdoit6 name a1/Re inc Ie H0+H1/L0 eps
			    (lambda (Ls/L0)
			      (lambda (result)
				(cond ((equal? result 'error) (display 1.0))
				      ((equal? result 'unstable) (display Ls/L0))
				      (else (display 0.0)))
				(display " "))))
		  (Ieloop (+ Ie dIe)))
		(begin (newline)
		       (iloop (+ inc di)))))
	  'done))))

(with-output-to-file "b18-10-80-20-80.dat"
  (lambda ()
    (HEMdoit-range6 18.0 0.3392 1.e-3 10. 80. 0.25 20. 80. 0.25)))
