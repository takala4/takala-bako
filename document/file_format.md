# ファイルフォーマット 

The english page is a below of japanese one.


## 行列データ

行列データをファイルに保持する場合，.csvファイル又は.datファイルを利用する．


### .csv

* 文字コード：UTF8
* 改行コード：CRLF
* 区切り文字：カンマ
* 原則としてBOM付き


```
a_{11}, ... , a_{1M}
  :   ,  .  ,  :
a_{N1}, ... , a_{NM}
```

### .dat

* 文字コード：UTF8
* 改行コード：CRLF
* 区切り文字：TAB


```
a_{11}     ...    a_{1M}
  :         .       :
a_{N1}     ...    a_{NM}
```


## ベクトルデータ

(数学的な意味の)ベクトルデータをファイルに保持する場合，.csvファイル又は.datファイルを利用する．


### .csv

```
b_{1}
:
b_{M}
```

### .dat

```
b_{1}
:
b_{M}
```


----

# File format 

Under constructed...


## Matrix Data


### .csv

.csv file is a delimited text file that uses a comma to separate values.

```
a_{11}, ... , a_{1M}
  :   ,  .  ,  :
a_{N1}, ... , a_{NM}
```

### .dat

.dat file is a delimited text file that uses a tab to separate values.

```
a_{11}, ... , a_{1M}
  :   ,  .  ,  :
a_{N1}, ... , a_{NM}
```

## Vector data


### .csv

```
b_{1}
:
b_{M}
```

### .dat

```
b_{1}
:
b_{M}
```
