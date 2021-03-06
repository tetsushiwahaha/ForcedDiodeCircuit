# 局所分岐集合追跡

## Files
* `bfFD.c` : 本体 
* `ns0.in` : 設定ファイル，Neimark Sacker 分岐の例
* `pd0.in` : 設定ファイル，周期倍分岐の例

## インストール

``` 
% make
```

## 実行

``` 
% ./bfFD ns0.in
```

## 設定ファイル内容

```
/* GAMMA, K, B0, B, X0, Y0 */
0.39820 0.62000 0.09000 0.29000 0.2925348 -1.0646094

/* パラメータ値，固定点 */

0.001 0.001 0.001 0.001		/* dk, dB0, dB */

/* 各パラメータの増分．符号を含めて指定．
下記の increment variable, variable for Newton Method の２つの
パラメータ空間において指定増分でスキャンする
*/

NS			
/* 分岐の種類． 
NS: Neimark-Sacker分岐，
PD or I: 周期倍分岐
G or T: 接線分岐
*/
/* 以下は FIX を参照 */
1e-12			/* eps */
100.0		/* explode */
1			/* number of period */
256			/* strip for runge-kutta */
20			/* iteration max */
1			/* extrapolation */
1.1			/* auto change parameter */
1			/* increment variable gamma = 0, k = 1, B0 = 2, B = 3 */
0	 		/* variable for Newton Method */
```

## 実行例
```
i% ./bfFD ns0.in 
i gamma  k    B0     B      x(1)     x(2)    mu1     mu2
3 0.39814 0.62000 0.09000 0.29000 0.2923012 -1.0641855 0.171+j0.985 (1.0000) o
4 0.39833 0.62100 0.09000 0.29000 0.2974157 -1.0669019 0.165+j0.986 (1.0000) o
3 0.39852 0.62200 0.09000 0.29000 0.3025256 -1.0696200 0.160+j0.987 (1.0000) o
3 0.39871 0.62300 0.09000 0.29000 0.3076310 -1.0723399 0.155+j0.988 (1.0000) o
3 0.39889 0.62400 0.09000 0.29000 0.3127321 -1.0750615 0.150+j0.989 (1.0000) o
3 0.39908 0.62500 0.09000 0.29000 0.3178288 -1.0777849 0.145+j0.989 (1.0000) o
```
Neimark-Sacker 分岐の例．
B0からdB0だけ徐々に変化させながら，gamma について解いた．
特性根の絶対値が1となっている．偏角も表示させるべきかもしれない．

```
% ./bfFD pd0.in 
i gamma  k    B0     B      x(1)     x(2)    mu1     mu2
3 0.38441 0.90200 0.09000 0.29000 1.8470323 -2.0209094 -1.00000 -1.00008 R
5 0.38399 0.90300 0.09000 0.29000 1.8489349 -2.0198568 -1.00000 -1.01381 R
3 0.38357 0.90400 0.09000 0.29000 1.8507913 -2.0187522 -1.00000 -1.02786 R
3 0.38314 0.90500 0.09000 0.29000 1.8525996 -2.0175938 -1.00000 -1.04223 R
3 0.38272 0.90600 0.09000 0.29000 1.8543579 -2.0163794 -1.00000 -1.05694 R
3 0.38229 0.90700 0.09000 0.29000 1.8560639 -2.0151065 -1.00000 -1.07201 R
3 0.38187 0.90800 0.09000 0.29000 1.8577153 -2.0137726 -1.00000 -1.08744 R
```
周期倍分岐の例．特性根の一つが -1 になっている．この例では他方の根もその
絶対値が1より大きいため，不安定な固定点であることがわかる．
