/*���������� ������� �ŷڱ��� : ���� 3-2 [ǥ 3.4] */

/* �ڷḦ �о�帮�� SAS ���� ����

����: SAS�� ��� ������ �����ݷ�(;)���� ������!!!!!!!!

(1) SAS dataset�� �̸� ����
    DATA �� ������  SAS dataset�� �̸��� �����Ѵ�. �Ʒ����� SAS dataset�� �̸��� liver�� �ȴ�.

(2) ������ ���� 
    INPUT �� �������� �о�帱 ������ �������  ����. ������������ �̸� �ڿ���  $�� ���δ�.
    ���� �� �ٿ� �� �� �̻��� ���������� ���� ���� �������� @@�� ���δ�.

(3) �ڷ�
    �����ʹ� ó���� CARDS; �� �����ϰ� INPUT ���� ������ ������ ������� �ڷḦ ��ġ�Ѵ�.
    ������ �����͵ڿ� �����ݷ�(;)�� ���δ�. 
*/

DATA test;
INPUT drug $ mi $ count;
CARDS;
1 1 73
1 2 18
2 1 141
2 2 196
;
PROC FREQ ;
TABLES drug*mi/MEASURES ;
WEIGHT count;
RUN;


/*ī������ ���� : ���� 3-3 [ǥ 3.6] */
DATA test1;
INPUT treat $ mi $ count;
CARDS;
a yes 139 
a no 10898
p yes 239
p no 10795
;
PROC FREQ DATA=test1;
WEIGHT count;
TABLE treat*mi/CHISQ ;
RUN;


/*�ƴϸ� ���� : ���� 3-5 [ǥ 3.14] */
DATA marriage;
INPUT before $ after $ count @@;
CARDS;
satisfy satisfy 23
satisfy unsatisfy 7
unsatisfy satisfy 18
unsatisfy unsatisfy 12
;

ODS SELECT McNemarsTest;
PROC FREQ order=data;
WEIGHT count;
TABLES before*after/AGREE;
RUN;


/*��ũ��-����-���� ���� : ���� 3-6 [ǥ 3.18] */
DATA hospital;
INPUT hospital $ trt $ recovery $ count @@;
CARDS;
A old yes 9 A old no   5
A new yes 11 A new no  6
B old yes 7 B old no   5
B new yes 8 B new no   3
C old yes 4 C old no   6
C new yes 7 C new no   5
D old yes 18 D old no  11
D new yes 26 D new no  4
;
PROC FREQ;
WEIGHT count;
TABLES hospital*trt*recovery/CMH NOROW NOCOL;
RUN;
