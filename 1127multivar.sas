proc import datafile = 'C:\Users\daess\desktop\income.csv'
dbms = csv replace out = income;
run;

proc print data = income;
run;

/*income 데이터에서 2번째 열 own을 분류*/
proc discrim data = income pool = yes out = local; /*pool = yes : 합동공분산을 쓰겠다.무조건 선형*/
class own;
var income size;
run;

proc discrim data = income pool = yes 
listerr out = local; /*pool = yes : 합동공분산을 쓰겠다.*/
class own;
var income size;
run;

proc print data = local;
run; /*결과창에서 _INTO_만 보면 됨.*/

data test; 
input xvalues own_test $ income size; 
cards; 
1 . 60.5 19.5 
2 . 79.0 18.0 
3 . 88.0 17.0 
;
proc print data = test;
run;

proc discrim data = local testdata = test testlist 
testout = test_out; 
class own; 
var income size; 
run; 

proc print data = test_out; 
run;

/*ex2*/
proc import datafile = 'C:\Users\daess\desktop\skull.csv'
dbms = csv replace out = skull;
run;
proc print data = skull (obs = 10);
run;

proc discrim data = skull pool = test listerr out = skull_out;
class group;
var x1--x4;
run;

proc discrim data = skull pool = test listerr out = skull_out simple;
class group;
var x1--x4;
priors prop; /*prior를 데이터 비율에 */
run;

/*cross validation*/
proc discrim data = skull pool = test  simple crosslisterr out = skull_out;
class group;
var x1--x4;
priors prop;
run;
