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
    in
        if f(str1, str2)
        then str1
        else str2
    end


val longest_string3 = longest_string_helper fn (x, y) => if x >= y then x else y

val longest_string4 = longest_string_helper fn (x, y) => if x > y then x else y


val longest_capitalized = fn strs => strs |> only_capitals |> longest_string1


val rev_string = fn str => str |> String.explode |> List.rev |> String.implode

fun first_answer f [] = raise NoAnswer
  | first_answer f x::xs =
    case f x of
	SOME v => v
      | NONE => first_answer f xs


      
fun all_answers f xs =
    let fun helper [] acc = SOME acc
          | helper x::xs acc =
                case f x of
                    NONE => NONE
                  | SOME v => helper xs (v @ acc)
    in helper xs []
    end


val count_wildcards = g (fn () => 1) (fn x => 0)

val count_wild_and_variable_lengths = g (fn () => 1) (fn str => String.size)

fun count_some_var (str pat) =
    g (fn () => 0) (fn str2 => if str == str2 then 1 else 0) pat


fun check_pat pat =
    let fun helper1 pat acc =
	    case pat of
		Variable x => x::acc
	      | TupleP ps => List.foldl(fn (p, acc) => helper p acc, acc, ps)
	      | _ => acc
	    
	fun helper2 [] = false
	  | helper2 x::xs = not List.exists (fn y => if x == y then true else false) xs andalso not helper2(xs)
												    
    in
	let val temp = helper1 pat []
	in 
	    case temp of 
		[] => true
	      | _ => helper2 temp
	end
    end


fun match (valu, pat) =
    case (valu, pat) of
	(_, Wildcard) => SOME []
     | (Unit, UnitP) => SOME []
     | (v, Variable s) => SOME [(s, v)]
     | (Const x, ConstP y) => if x = y then SOME [] else NONE
     | (Tuple vs, TupleP ps) => if List.length vs <> List.length ps
			      then NONE
			      else 
				  let pairList = List.pair(vs, ps)
				  in 
				      case all_answers (fn (x, y) => match(x, y)) pairList of
					| NONE => NONE
					| SOME lst => SOME lst
				  end
     | (Constructor(s1, v), Constructor(s2, p)) => if s1 = s2
						   then match(v, p)
						   else NONE
     | (_, _) => NONE

