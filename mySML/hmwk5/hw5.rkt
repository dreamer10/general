;; Programming Languages, Homework 5

#lang racket

(require racket/trace)

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
  (if (aunit? mupl_list)
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
        
        [(int? e) e]
        
        [(closure? e) e]
        
        [(aunit? e) e]
        
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
         (closure env e)]

        [(call? e)
         (let ([fun_closure (eval-under-env (call-funexp e) env)]
               [fun_arg (eval-under-env (call-actual e) env)])
           (if (not (closure? fun_closure))
               (error "MUPL no valid function in the environment")
               (letrec ([c_env (closure-env fun_closure)]
                     [func (closure-fun fun_closure)]
                     [func_name (fun-nameopt func)]
                     [func_arg (fun-formal func)]
                     [func_body (fun-body func)])
                 (if func_name
                     (eval-under-env func_body
                                     (cons (cons func_name fun_closure) (cons (cons func_arg fun_arg) c_env)))
                     ; else if it's anonymous function
                     (eval-under-env func_body (cons (cons func_arg  fun_arg) c_env))))))]
        
        [(mlet? e)
         (let ([e1_val (eval-under-env (mlet-e e) env)])
           (eval-under-env (mlet-body e) (cons (cons (mlet-var e) e1_val) env)))]
        
        [(apair? e)
         (let ([v1 (eval-under-env (apair-e1 e) env)]
               [v2 (eval-under-env (apair-e2 e) env)])
           (apair v1 v2))]
        
        [(fst? e)
         (let ([v (eval-under-env (fst-e e) env)])
           (if (apair? v)
               (apair-e1 v)
               (error "MUPL fst applied not to a pair")))]
        [(snd? e)
         (let ([v (eval-under-env (snd-e e) env)])
           (if (not (apair? v))
               (error "MUPL snd applied NOT to a pair")
               (apair-e2 v)))]
        
        [(isaunit? e)
         (let ([v (eval-under-env (isaunit-e e) env)])
           (if (aunit? v)
               (int 1)
               (int 0)))]
                                 
                                                      
        [#t (error (format "bad MUPL expression: ~v" e))]))

;; Do NOT change
(define (eval-exp e)
  (eval-under-env e null))
        
;; Problem 3

(define (ifaunit e1 e2 e3) (ifgreater (isaunit e1) (int 0) e2 e3))

(define (mlet* lstlst e2)
  (if (null? lstlst)
      e2
      (let ([pr (car lstlst)])
        (mlet (car pr) (cdr pr) (mlet* (cdr lstlst) e2)))))


; Other's code as reference
(define (ifeq e1 e2 e3 e4)
  (mlet* (list (cons "_x" e1) (cons "_y" e2))
         (ifgreater (var "_x") (var "_y") e4
                    (ifgreater (var "_y") (var "_x") e4
                               e3))))


#| My wrong code:

(define (ifeq e1 e2 e3 e4)
  (let* ([pr (apair e1 e2)]
         [hd (fst pr)]
         [tl (snd pr)]
         [thunk_e3 (fun #f "_x" e3)]
         [thunk_e4 (fun #f "_y" e4)])
     (mlet* (list (cons "_x" (call thunk_e3 (int 0))) (cons "_y" (call thunk_e4 (int 0))))
         (ifgreater (var "_x") (var "_y") e4
                    (ifgreater (var "_y") (var "_x") e4 e3)))))

|#



(define mupl-map
  (fun #f "f" 
       (fun "my_map" "lst"
       (ifeq (isaunit (var "lst")) (int 1)
             (aunit)
             (apair (call (var "f") (fst (var "lst"))) (call (var "my_map") (snd (var "lst"))))))))


; When debugging, this is ##VERY IMPORTANT##:
(trace eval-under-env)


#|
Other's code as reference:

; Problem 4
(define mupl-map
  (fun #f "func"
       (fun "iter" "lst"
            (ifaunit (var "lst") (aunit)
                     (apair (call (var "func") (fst (var "lst")))
                            (call (var "iter") (snd (var "lst"))))))))
|#


(define mupl-mapAddN 
  (mlet "map" mupl-map
        (fun #f "i"
             (call (var "map") (fun #f "x" (add (var "x") (var "i")))))))

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
