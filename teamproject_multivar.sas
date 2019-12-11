/*주성분분석*/
proc import datafile = 'C:\Users\daess\Desktop\newdata.xls'
dbms = xls out = c replace;
run;

proc print data = c;
run;

proc standard data = c mean = 0 std = 1 out = bb;
var y x1--x41;
run;

proc print data = bb;
run;

proc princomp data = bb out = bbb;
var x41 x40 x5 x6 x8 x37 x38;
run;

/* 데이터 내보내기 */

proc export data = bbb
outfile = 'C:\Users\daess\Desktop\prin.csv'
dbms = csv replace;
run;


/* ÆÇº°ºÐ¼® */

proc import datafile =  'C:\Users\daess\Desktop\realprin.csv'
dbms = csv replace out = test;
run;

proc discrim data = test listerr pool = test manova out = test_out;
class group4;
var prin1 prin2 prin3;
priors prop;
run;

proc import datafile =  'C:\Users\daess\Desktop\realprin.csv'
dbms = csv replace out = test;
run;

proc discrim data = test listerr pool = test manova out = test_out;
class group5;
var prin1 prin2 prin3;
priors prop;
run;

proc import datafile =  'C:\Users\daess\Desktop\realprin.csv'
dbms = csv replace out = test;
run;

proc discrim data = test listerr pool = test manova out = test_out;
class group6;
var prin1 prin2 prin3;
priors prop;
run;

proc import datafile =  'C:\Users\daess\Desktop\realprin.csv'
dbms = csv replace out = test;
run;

proc discrim data = test listerr pool = test manova out = test_out;
class group7;
var prin1 prin2 prin3;
priors prop;
run;


/* 데이터 */

proc import datafile =  'C:\Users\daess\Desktop\prin.csv'
dbms = csv replace out = aaa;
run;

proc print data = aaa;
run;

proc standard data = aaa mean = 0 std = 1 out = aa_;
var aaa;
run;

proc print data = aa_;
run;

proc export data = aa_ outfile = 'C:\Users\daess\desktop\newnewprin.csv' dbms = csv
replace;
run;
