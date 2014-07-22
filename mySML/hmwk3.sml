(* Coursera Programming Languages, Homework 3, Provided Code *)

exception NoAnswer

datatype pattern = Wildcard
		 | Variable of string
		 | UnitP
		 | ConstP of int
		 | TupleP of pattern list
		 | ConstructorP of string * pattern

datatype valu = Const of int
	      | Unit
	      | Tuple of valu list
	      | Constructor of string * valu

fun g f1 f2 p =
    let 
	val r = g f1 f2 
    in
	case p of
	    Wildcard          => f1 ()
	  | Variable x        => f2 x
	  | TupleP ps         => List.foldl (fn (p,i) => (r p) + i) 0 ps
	  | ConstructorP(_,p) => r p
	  | _                 => 0
    end

(**** for the challenge problem only ****)

datatype typ = Anything
	     | UnitT
	     | IntT
	     | TupleT of typ list
	     | Datatype of string

(**** you can put all your code here ****)


fun only_capitals [] = []
  | only_capitals strs = 
    List.filter(fn str => Char.isUpper(String.sub str), strs)

fun longest_string1 [] = ""
  | longest_string1 strs =
    List.foldl(fn (str1, str2) => if String.size str2 > String.size str1
				  then str2
				  else str1,
	       "", strs)

fun longest_string2 [] = ""
  | longest_string2 strs = 
    List.foldl(fn (str1, str2) => if String.size str2 >= String.size str1
				  then str2
				  else str1,
	       "",
	       strs)


fun longest_string_helper f [] = ""
  | longest_string_helper f str1::xs =
    let val str2 = longest_string_helper f xs
    if f(str1, str2)
    then str1
    else str2


val longest_string3 = longest_string_helper fn (x, y) => if x >= y then x else y

val longest_string4 = longest_string_helper fn (x, y) => if x > y then x else y

