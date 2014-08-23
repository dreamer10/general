;; Programming Languages, Homework 5

#lang racket
(provide (all-defined-out)) ;; so we can put tests in a second file

;; definition of structures for MUPL programs - Do NOT change
(struct var  (string) #:transparent)  ;; a variable, e.g., (var "foo")
(struct int  (num)    #:transparent)  ;; a constant number, e.g., (int 17)
(struct add  (e1 e2)  #:transparent)  ;; add two expressions
(struct ifgreater (e1 e2 e3 e4)    #:transparent) ;; if e1 > e2 then e3 else e4
(struct fun  (nameopt formal body) #:transparent) ;; a recursive(?) 1-argument function
(struct call (funexp actual)       #:transparent) ;; function call
(struct mlet (var e body) #:transparent) ;; a local binding (let var = e in body) 
(struct apair (e1 e2)     #:transparent) ;; make a new pair
(struct fst  (e)    #:transparent) ;; get first part of a pair
(struct snd  (e)    #:transparent) ;; get second part of a pair
(struct aunit ()    #:transparent) ;; unit value -- good for ending a list
(struct isaunit (e) #:transparent) ;; evaluate to 1 if e is unit else 0

;; a closure is not in "source" programs; it is what functions evaluate to
(struct closure (env fun) #:transparent) 

;; Problem 1

(define (racketlist->mupllist rkt_list)
  (if (null? rkt_list)
      (aunit)
      (apair (car rkt_list) (racketlist->mupllist (cdr rkt_list)))))

;; Problem 2

(define (mupllist->racketlist mupl_list)
  (if (isaunit mupl_list)
      null
      (cons (apair-e1 mupl_list)  (mupllist->racketlist (apair-e2 mupl_list)))))

;; lookup a variable in an environment
;; Do NOT change this function
(define (envlookup env str)
  (cond [(null? env) (error "unbound variable during evaluation" str)]
        [(equal? (car (car env)) str) (cdr (car env))]
        [#t (envlookup (cdr env) str)]))

;; Do NOT change the two cases given to you.  
;; DO add more cases for other kinds of MUPL expressions.
;; We will test eval-under-env by calling it directly even though
;; "in real life" it would be a helper function of eval-exp.
(define (eval-under-env e env)
  (cond [(var? e)
         (envlookup env (var-string e))]
        [(add? e) 
         (let ([v1 (eval-under-env (add-e1 e) env)]
               [v2 (eval-under-env (add-e2 e) env)])
           (if (and (int? v1)
                    (int? v2))
               (int (+ (int-num v1) 
                       (int-num v2)))
               (error "MUPL addition applied to non-number")))]
        ;; CHANGE add more cases here
        [(ifgreater? e)
         (let ([v1 (eval-under-env (ifgreater-e1 e) env)]
               [v2 (eval-under-env (ifgreater-e2 e) env)])
               (if (and (int? v1)
                        (int? v2))
                   (if (> (int-num v1) (int-num v2))
                       (eval-under-env (ifgreater-e3 e) env)
                       (eval-under-env (ifgreater-e4 e) env))
                   (error "MUPL comparison applied to non-number")))]
        [(fun? e)
         (let ([fun_name (fun-nameopt e)])
           (if fun_name
               (closure (cons (cons fun_name (fun-body e)) env) (cons (fun-nameopt e) (fun-formal e)))
               (closure (cons (cons (fun-formal e) (fun-body e)) env) (cons (fun-nameopt e) (fun-formal e)))))]
        
        [(call? e)
         (let ([fun_closure (eval-under-env (call-funexp e) env)]
               [fun_arg (eval-under-env (call-actual e))])
           (if (not (closure? fun_closure))
               (error "MUPL no valid function in the environment")
               (let* ([c_env (car fun_closure)]
                     [fun_pr (cdr fun_closure)]
                     [fun_name (car fun_pr)]
                     [fun_argname (cdr fun_pr)])
                 (if fun_name
                     (eval-under-env (envlookup c_env fun_name)
                                     (cons (fun_name fun_closure) (cons (fun_argname fun_arg) c_env)))
                     ; else if it's anonymous function
                     (eval-under-env (envlookup c_env fun_argname) c_env)))))]
        
        [(mlet? e)
         (let ([e1_val (eval-under-env (mlet-e e) env)])
           (eval-under-env (mlet-body e) (cons ((mlet-var e) e1_val) env)))]
        
        [(apair? e)
         (let ([v1 (eval-under-env (apair-e1 e) env)]
               [v2 (eval-under-env (apair-e2 e) env)])
           (apair v1 v2))]
        
        [(fst? e)
         (let ([v (eval-under-env e env)])
           (if (apair? v)
               (apair-e1 e)
               (error "MUPL fst applied not to a pair")))]
        [(snd? e)
         (let ([v (eval-under-env e env)])
           (if (not (apair? v))
               (error "MUPL snd not applied to a pair")
               (apair-e2 v)))]
        
        [(isaunit? e)
         (let ([v (eval-under-env e)])
           (if (aunit? v)
               (int 1)
               (int 0)))]
                                 
                                                      
        [#t (error (format "bad MUPL expression: ~v" e))]))

;; Do NOT change
(define (eval-exp e)
  (eval-under-env e null))
        
;; Problem 3

(define (ifaunit e1 e2 e3) (ifgreater (isaunit e1) 0 e2 e3))

(define (mlet* lstlst e2)
  (if (null? lstlst)
      e2
      (let ([pr (car lstlst)])
        (mlet (car pr) (cdr pr) (mlet* (cdr lstlst) e2)))))

  
(define (ifeq e1 e2 e3 e4)
  (let* ([pr (apair e1 e2)]
         [hd (fst pr)]
         [tl (snd pr)]
         [thunk_e3 (fun #f "x" e3)]
         [thunk_e4 (fun #f "y" e4)])
    (if (and (int? hd) (int? tl))
        (if (equal? (int-num hd) (int-num tl))
            (call thunk_e3 0)
            (call thunk_e4 0))
        (error "equal? applied to non-intergers"))))
            
    
;; Problem 4

(define mupl-map "CHANGE")

(define mupl-mapAddN 
  (mlet "map" mupl-map
        "CHANGE (notice map is now in MUPL scope)"))

;; Challenge Problem

(struct fun-challenge (nameopt formal body freevars) #:transparent) ;; a recursive(?) 1-argument function

;; We will test this function directly, so it must do
;; as described in the assignment
(define (compute-free-vars e) "CHANGE")

;; Do NOT share code with eval-under-env because that will make
;; auto-grading and peer assessment more difficult, so
;; copy most of your interpreter here and make minor changes
(define (eval-under-env-c e env) "CHANGE")

;; Do NOT change this
(define (eval-exp-c e)
  (eval-under-env-c (compute-free-vars e) null))
