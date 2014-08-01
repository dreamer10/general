fun is_older (date1 : int * int * int, date2 : int * int * int) = 
    if (#1 date1) = (#1 date2) andalso (#2 date1) = (#2 date2) andalso (#3 date1) < (#3 date2)
    then true
    else if (#1 date1) = (#1 date2) andalso (#2 date1) < (#2 date2)
    then true
    else if (#1 date1) < (#1 date2)
    then true
    else false


fun number_in_month (dates : (int * int * int) list, month : int) = 
    if null dates
    then 0
    else if (#2 (hd dates)) = month
    then 1 + number_in_month(tl dates, month)	    
    else number_in_month(tl dates, month)


fun number_in_months (dates : (int * int * int) list, months : int list) =
    if null months
    then 0
    else number_in_month(dates, (hd months)) + number_in_months(dates, (tl months))


fun dates_in_month (dates : (int * int * int) list, month : int) = 
    if null dates
    then []
    else if (#2 (hd dates)) = month
    then (hd dates) :: dates_in_month(tl dates, month)
    else dates_in_month(tl dates, month)

fun dates_in_months (dates : (int * int * int) list, months : int list) = 
    if null months
    then []
    else dates_in_month(dates, hd months) :: dates_in_months(dates, tl months)





(* fun dates_in_months (dates : (int * int * int) list, months : int list) = *)
(*     let fun isInMonthList(monthOfDate : int, months: int list) =  *)
(* 	    if null months *)
(* 	    then false *)
(* 	    else if monthOfDate = #2 (hd months) *)
(* 	    then true *)
(* 	    else isInMonthList(monthOfDate, tl months) *)
(*     in let fun getList(dates : (int * int * int) list, months : int list) =  *)
(* 	       if null dates *)
(* 	       then [] *)
(* 	       else if isInMonthList((#2 (hd dates)), months) *)
(* 	       then (hd dates) :: getList(tl dates, months) *)
(* 	       else getList(tl dates, months) *)
(*        in getList(dates, months) *)
(*        end *)
(*     end *)

fun get_nth (strList : string list, n : int) = 
    if n = 1
    then hd strList
    else get_nth(tl strList, n - 1)

fun date_to_string (date : int * int * int) = 
    let val monthList = ["January", "Februray", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    in Int.toString(#1 date) ^ get_nth(monthList, #2 date) ^ Int.toString(#3 date)
    end


fun number_before_reaching_sum (sum : int, posNums : int list) =
    let fun helper (sum : int, posNums : int list, ret : int) =
	    if sum < (hd posNums)
	    then ret
	    else helper(sum - (hd posNums), tl posNums, ret + 1)
    in helper (sum, posNums, 0)
    end

fun what_month (day : int) =
    let val days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    in number_before_reaching_sum(day, days_in_month) + 1
    end

fun month_range (day1 : int, day2 : int) = 
    if day1 > day2
    then []
    else 
	let fun count_from_to (mon1 : int, mon2 : int) =
		if mon1 = mon2
		then mon1 :: []
		else mon1 :: count_from_to(mon1 + 1, mon2)
	in count_from_to(what_month(day1), what_month(day2))
	end



fun oldest (dateList : (int * int * int) list) =
    if null dateList
    then NONE
    else
	let fun getOldestDate (dateList : (int * int * int) list) =
		if null (tl dateList)
		then hd dateList
		else
		    let val tmp = getOldestDate(tl dateList)
		    in if is_older(hd dateList, tmp)
		       then hd dateList
		       else tmp
		    end
	in 
	    SOME (getOldestDate dateList)
	end
