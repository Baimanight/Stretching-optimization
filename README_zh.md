# Stretching-optimization


[English](README.md) | [中文](README_zh.md)

##### 本文针对[MSNoise1.6](http://www.msnoise.org) 中stretching代码部分进行讨论，以及考虑钟差校正对ST算法进行优化。目前只是一个比较简单的尝试，代码部分较为粗糙，一些想法也有待验证。但是只要路径正确，就应该不会报错。。。

### 使用方法

ST_FPC（Fold-Prediction-Correction），可用于替代MSNoise流程中Stack后的stretching计算部分。与[MSNoise](http://www.msnoise.org)结合使用，利用STACKS做为输入，结果输出到Str文件夹，参数需要在ST_FPC.m代码头部设置（文件路径、流逝时间窗口、拉伸范围、步长等）。

* readsac.m 为SAC官方函数，需要与ST_FPC.m放于同一文件夹下。

* stretch.py 为MSNoise计算ST函数。

* MATLAB需要保证python环境，以调用map_coordinates函数。

### 算法解释

个人认为，MWCS（移动窗互谱法）在原理上是明显优于ST（压缩拉伸法）的。在计算参考波形和叠加后波形的dt/t，即两者拉伸系数时，所最为依赖的是波形间的相位关系。MWCS是直接依靠相位关系来拟合得到dt/t，而ST则是通过可以反应相位相似程度的CC（相关系数）来认定dt/t(Obermann and Hillers, 2019)。后者经历了由波形到CC、CC到dt/t的两层失真，对原波形信息的提取是不如MWCS的。但实际计算中，ST涉及参数较少，使用较为方便，或可作为参考。

由经验格林函数计算速度变化 $dt/t$ ，其实就是计算Ref与Days波形拉伸系数，假设两波形分别对应函数 $f（x）$与 $f（ax+b）$， $a$就是拉伸系数，而“钟差”在波形中的体现就是 $f（x）$与 $f（ax+b）$中的 $b$.  故将分别从dt/t和offset两节介绍，测试数据为某六个相近台站多月的连续波形数据。

### dt/t

将Ref波形按照设置的拉伸区间依步长变化拉伸，计算所有拉伸Ref波形与Days波形的互相关系数Correlation Coefficient（Pearson），取CC最大时对应的拉伸系数做为dt/t。

图一是一段互相关函数波形-120\~120s，黄色部分为选取的流逝时间窗口20\~80s。拉伸时，原代码选用黄色部分数据点（包括置零的非窗口区），拉伸时间横轴再利用scipy.ndimage.map_coordinates 函数进行插值。计算CC时，使用拉伸后-120~120s全部数据点与非窗口区置零的Days数据点计算CC。

![](Figure/c755c4a0-e22c-11ef-b911-b3d360b4824a.jpeg?v=1&type=image)
（图一）


首先对原python原代码进行复现，MATLAB中的插值函数interp1相较于map_coordinates缺少预滤波过程，效果不好，所以直接调用了python中的该函数。结果如图二所示，颜色代表CC大小（余图一致），基本复现原python代码。

![](Figure/8dbdb1c0-e2d7-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image)
（图二）


 原代码似乎存在以下问题：

当压缩波形时，窗口边缘数值需要进行外插而默认为0，但已知全段波形，完全可以加长区段以保证区段外插全为有效值；

而当拉伸波形时，区段外会产生额外的数据点，这部分数据点原代码并没有清除，而是直接参与CC计算，并且两段数据间存在置零区间，对应的测试部分为[-20, 20]，这些都将影响最终CC结果。



对于波形区段，我选择将波形对折取平均，选取窗口内的数据。结果如图三，仅通过对区段的选择，便使89%的CC有了平均约10%的提高。

![](Figure/31db0ce0-e2d7-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image)
（图三）


对于外插问题，可以通过简单计算加长区段，即可保证外插部分全为有效值 。图四所示，Fold Width为加长区段的结果，可见85%CC均有小幅提升。

![](Figure/0ce2d0f0-e235-11ef-b911-b3d360b4824a.jpeg?v=1&type=image)
（图四）




### Offset

不同于MWCS(Stehly et al. , 2007)，ST并没有可靠的方式直接处理钟差，只能同样由CC大小来认定可能的钟差，我将这样得到的“钟差”称之为offset（类钟差的补偿）。一个简单的想法是，通过分别拟合得到Ref与Days的函数表达式，再对比函数系数找到 f（x） 与f（x+b）的b。如果直接聚焦于选定的流逝窗口，这将是一件十分困难的事。但是通过观察整体波形特征，并进行了一些试算验证，摸索出了这种基于Gauss拟合的预报—校正方法。   &#x20;

如图，黄色和紫色波形分别为同一台站对Day和Ref的互相关函数，对数据取绝对值后，选取四阶高斯拟合，得到的系数以加权平均的方式计算中心，将两波形拟合的中心之差作为offset。


![](Figure/01d53ce0-e2d8-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image)
（图五）


图六中，利用Gauss4拟合所得到的offset呈现出一定横向分布特性，也就是与台站对的相关性，而钟差本身也就与各台站对有着对应关系，这样得到的offset或可认定为“钟差”

![](Figure/9a4db2d0-e473-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image)
（图六）


图七中纵轴为15个台站对，可见横向颜色相近条带，指示offset结果与台站对呈一定相关性

![](Figure/0b0b7de0-e46f-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image)
（图7）


##### 但此时的offset与波形相位关系并无直接联系，只是模拟的结果，且仅依据此时offset得到的CC结果仅半数略微改善，并不可靠，还需要结合CC来选取更为合适的offset

![](Figure/1f9fcc70-e52d-11ef-baf1-a700618c26d1.jpeg?v=1&type=image)
（图八）


此外，Gauss拟合时依赖适当的参数选取，不当的参数对结果影响较大。如图，当不设置拟合初始值时，offset分布十分不均，存在明显误差。

![](Figure/438147e0-e474-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image)
（图九）


##### 可如果只依靠CC的最大值作为offset，在设定范围内offset又会呈现出随机分布的特征，如图十。也就是说，无法完全依赖CC作为有效的参考指标。

![](Figure/1c80f5d0-e886-11ef-8864-3d63cc015e26.jpeg?v=1&type=image)
（图十）


再次观察其与台站分布相关性（图十一），注意到大量颜色交错带，即数据剧烈变化的振荡区域，此段offset一定不可认为钟差。

![](Figure/83dde800-e52a-11ef-baf1-a700618c26d1.jpeg?v=1&type=image)
（图十一）


受此启发，可识别出振荡区域，并将该区内不可靠的offset全部置零处理。加之边界限制，可见offset有了不小的改善。

![](Figure/96fe9370-e887-11ef-8864-3d63cc015e26.jpeg?v=1&type=image)
（图十二）


#### 针对上述问题，引入了预报—校正系统，流程如图十三

#### 1. 预报：先由Gauss4 拟合预报出offset大致分布，由此预报拉伸系数pre-dt/t

#### 2. 校正：依据预报拉伸系数拉伸波形，在预报offset范围内依据CC最大值得到校正的offset，剔除offset异常区域

#### 3. 最后利用校正的offset微调窗口，计算CC、取其最大时对应的拉伸系数作为Corr-dt/t

![](Figure/624f6590-eb7b-11ef-bb64-49bc940d885a.jpeg?v=1&type=image)
（图十三）


##### 至此考虑“钟差”校正的优化结果如图，相较Fold的结果，88%的CC有了一定提高

![](Figure/5ae41860-e2c4-11ef-b8f7-7fbcf1303d2b.jpeg?v=1&type=image)
（图十四）


##### 最终，相较于原代码，97%的CC结果都有了平均约15%的提高。

![](Figure/cf4b1880-e2c3-11ef-b8f7-7fbcf1303d2b.jpeg?v=1&type=image)
（图十五）


##### 而CC的提高，让依据CC的加权平均成为可能， 所示仅为三点的加权平均结果，但速度变化已经十分稳定

![](Figure/914b5820-eba5-11ef-b0c0-5937d7efb862.jpeg?v=1&type=image)
（图十六）





### 问题总结

##### 参考MSNoise1.6的代码， 通过Fold改变参与计算的波形，以及Gauss拟合、预报—校正系统和对振荡数据的筛除，最终在一定程度上改进了ST算法，提高了dt/t的CC。&#x20;

##### 但是仍然存在如下一些问题，&#x20;

* ST的offset与mwcs的“钟差”相比如何？

* 如何有效地选取合适的Gauss拟合参数或有什么更好的拟合方法？

* offset振荡区除了置零是否还有更好的方法？

* offset的分布可见一定纵向相关，或含有有某种周期性，是否可以拟合用以校正？

* CC的选择上Pearson系数似乎并不能很好适用此问题，距离相关系数或MIC的效果如何呢？         &#x20;

*

##### 这些都需要进一步的研究，个人学习背景噪声时间较短，对相关问题理解尚不深入，对于本文中存在的任何问题，期待各位指正！

***

References

[1]Obermann A, Hillers G. Chapter Two-Seismic time-lapse interferometry across scales [J]. Advance in Geophysics, 2019, 60: 65-143.

[2]Stehly L., Campillo M., Shapiro N. M. 2007. Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. Geophys. J. Int., 171(4): 223-230.&#x20;
