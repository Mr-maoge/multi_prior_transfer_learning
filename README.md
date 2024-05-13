# 论文

基于多源域知识迁移学习的小微企业信用评分

[基于多源域知识迁移学习的小微企业信用评分 - 中国知网 (cnki.net)](https://kns.cnki.net/kcms2/article/abstract?v=WdAl4K16JyWnrwJ7JhdrFG9Hm0-CbGTRZvjGKjefiabQpA5_3qK_5xUSWYVNjvWGtADvd9DZH0MQXUGJhasHLaN29iHMq0f9coM4hJpktAt6rTiMuJWqxL9eA9Oq9J8wp25IytxvTMGCaVLrjjQN3A==&uniplatform=NZKPT)



# 摘要

针对新业务, 新场景下金融机构目标数据集 “高维小样本” 的问题, 本文提出了基于 多源域知识迁移学习的小微企业信用风险测度方法, 其能够迁移学习其它数据源 (源域) 的 知识以提升目标域模型的预测效果. 该方法通过对来自多个源域的多种源域知识进行归纳提 取, 进而将其纳入目标域模型的构建中, 可以充分利用源域知识, 提升目标域模型的估计精 度. 另外, 模型无需获取各源域的原始数据, 因此很大程度上降低了数据传输中隐私泄露的 风险. 模拟实验和企业信用评分的实例数据验证了所提方法的可行性及其在变量选择, 系数 估计和分类预测上的良好效果. 该方法能够在隐私限制的背景下有效迁移源域知识以克服信 用评分中目标数据集信息量不足, 而导致估计效果较差的问题.



# 联系人

李晶茂，[mr_ljm@foxmail.com](mailto:mr_ljm@foxmail.com)



# 代码文件主要函数和运行方法说明

## 主要文件

* dgp.R

  * data_generator1

    用于生成模拟数据

    

* Method.R

  * logistic.CD

    使用坐标下降法求解惩罚Logistic回归

  * cv.logistic

    基于logistic.CD，使用交叉验证法挑选调节参数

  * transfer.logistic

     多源域知识迁移学习方法（本文方法）

  * cv.transfer.logistic

    基于transfer.logistic，使用交叉验证法挑选调节参数

    

* simulation_funcs.R

  * simulation1

    用于模拟实验中情形 1和情形 5的结果计算函数

    （各源域知识均为系数估计值的场景）

  * simulation2

    用于模拟实验中情形 2的结果计算函数

    （各源域知识均为显著变量集的场景）

  * simulation3

    用于模拟实验中情形 3的结果计算函数

    （各源域知识包含系数估计值和显著变量集两种形式的场景）

  * simulation4

    用于模拟实验中情形 3的结果计算函数

    （各源域知识的质量均较差的场景）

    

* run_simulation_examples.R

  使用并行计算计算各个不同模拟情形、参数设定下100次模拟的结果。



* other_setting3.csv

  存储了模拟实验中参数的不同组合设定，具体包括：

  | rou  | n    | p    |
  | ---- | ---- | ---- |
  | 0.3  | 80   | 60   |
  | 0.5  | 80   | 60   |
  | 0.3  | 120  | 80   |
  | 0.5  | 120  | 80   |
  | 0.3  | 160  | 100  |
  | 0.5  | 160  | 100  |



* /result/

  结果文件夹

  

## 运行方法

运行 `run_simulation_examples.R` 文件，所以模拟结果将存储于`/reuslt` 文件夹中。