# FastNSA Python 中文使用指南

FastNSA 是一个由 C++ 后端加速的网络约束空间点分析工具包，当前 Python
接口主要提供以下能力：

- 将点事件匹配到道路边，得到 edge_id 和沿边 offset；
- 网络核密度估计（NKDE/NKDV）；
- Network K-function；
- Spatiotemporal Network K-function；
- Monte Carlo 理论参考、包络和模拟 p-value；
- CSV、DataFrame、NPZ 与 xarray 结果输出。

本文面向 FastNSA 0.3.x。更完整的参数签名可查阅
[英文 API reference](api-reference.md)。

## 安装

FastNSA 要求 Python 3.11 或更高版本，基础包只依赖 NumPy：

~~~bash
python -m pip install fastnsa
python -m fastnsa._selftest
~~~

具体可用平台和 Python 版本由每次发布的 wheel 文件决定。如果找不到匹配
的 wheel，pip 可能尝试源码编译，此时需要 CMake 和支持 C++17 的编译器。

按实际需求安装可选依赖：

~~~bash
python -m pip install "fastnsa[dataframe]"  # pandas 输出
python -m pip install "fastnsa[xarray]"     # ST 二维结果
python -m pip install "fastnsa[geo]"        # GeoPandas 与 Shapely
python -m pip install "fastnsa[osm]"        # NetworkX 与 OSMnx
python -m pip install "fastnsa[all]"        # 全部可选依赖
~~~

如果在 QGIS 中使用，必须安装到 QGIS 实际使用的 Python 环境。系统 Python、
Conda 环境和 QGIS Python 可能是三个相互独立的环境。

## 最小可运行示例

下列示例不依赖 GeoPandas、NetworkX 或 OSMnx，可用于检查安装是否正确：

~~~python
import fastnsa as nsa

with nsa.make_line_network(length=1000.0) as network:
    with nsa.make_line_events(
        network,
        offsets=[100.0, 350.0, 700.0, 900.0],
    ) as events:
        result = nsa.network_k_function(
            network=network,
            points=events,
            r_values=[50.0, 100.0, 250.0],
            method=nsa.KFunctionMethod.STARS,
            simulations=5,
            random_state=2026,
        )
        print(result.summary())
        result.to_csv("network_k.csv")
        result.save_npz("network_k.npz")
~~~

5 次模拟适合功能测试，不适合正式统计推断。正式分析应根据所需 p-value
分辨率增加模拟次数，并记录随机种子。

## 理解核心数据结构

FastNSA 的分析函数都接收一个网络和一组已经匹配到网络上的点。

| 低依赖接口 | 高层接口 | 用途 |
| --- | --- | --- |
| NativeNetwork | Network | 保存道路拓扑、边长度和边几何 |
| NativePointSet | PointEvents | 保存 edge_id、offset 和可选时间 |

NativeNetwork 和 NativePointSet 只需要 NumPy 与原生库，适合 QGIS 插件或
其它宿主程序直接传递数组。Network 和 PointEvents 提供 NetworkX、
GeoPandas、CRS 和 OSMnx 兼容功能。

原生对象持有 C++ 内存，推荐使用 with，或者按照“先关闭 events，再关闭
network”的顺序显式调用 close。

## CRS 与距离单位

这是使用 FastNSA 时最重要的正确性要求：

- 道路网络必须使用适合研究区域的投影 CRS；
- 道路边长、bandwidth、lixel length、r 和 s 必须使用相同单位；
- FastNSA 不会自动把经纬度的“度”转换成地表米制距离；
- 重新投影、修改道路或改变边顺序后，旧 edge_id/offset 必须作废并重新匹配。

例如使用米制投影道路时，bandwidth=500 表示 500 米；直接使用 EPSG:4326
时则会被当作 500 个坐标单位，而不是 500 米。

## 加载真实点数据

安装 geo 可选依赖后，可以加载 CSV、GeoPackage、GeoJSON 和 GeoParquet：

~~~python
points = nsa.load_points(
    "events.gpkg",
    layer="events",
    time_col="event_time",
)
~~~

CSV 默认读取 lon 和 lat 字段，并赋予 EPSG:4326：

~~~python
points = nsa.load_points(
    "events.csv",
    x_col="longitude",
    y_col="latitude",
    crs="EPSG:4326",
)
~~~

输入必须是非空 Point 图层且包含 CRS。加载后仍应把点匹配到已经投影的网络，
PointEvents 会在点和网络 CRS 不一致时先将点转换到网络 CRS。

## 构建网络和匹配事件

点到道路的最近邻匹配方法来源于 Ye 等人在 KDD 2026 发表的 *An Efficient
and Accurate Grid Compression Solution for Point-to-Line Nearest Neighbor
Search*。论文与实现链接见[算法参考文献](references.md#point-to-line-nearest-edge-matching)。

通过 OSMnx 构建网络：

~~~python
network = nsa.Network.from_osm(
    place="Gainesville, Florida, USA",
    network_type="drive",
)
events = nsa.PointEvents(
    points,
    target_network=network,
    mapping_method=nsa.MappingMethod.RTREE,
    n_threads=4,
)
~~~

R-tree 是一般用途的默认匹配方法。comp_e 可作为小数据集上的穷举对照，
comp_a 使用 comp_eps 容差。真实项目应在目标道路规模上比较精度和耗时。

## 复用 nearest-edge 预处理结果

如果点图层已经包含 edge_id 和 offset，并且可以确认道路网络没有变化，后续
分析无需再次匹配：

~~~python
with nsa.preprocess_events(
    points,
    network,
    edge_id_col="edge_id",
    offset_col="offset",
    time_col="event_time",
) as prepared:
    result = nsa.network_k_function(
        network=network,
        points=prepared.events,
        r_max=5000.0,
        r_step=100.0,
    )
    print(prepared.reused_mapping)
    print(prepared.timings)
~~~

建议同时保存 network signature。仅比较 CRS 或图层名称不足以证明映射仍然
有效，边几何、拓扑和内部顺序都必须一致。

## Network KDE

ADA/IA 增广方法和 LION 分别来源于 PVLDB 2021 与 PVLDB 2024 论文，详见
[NKDE 参考文献](references.md#network-kernel-density-visualization)。

~~~python
kde = nsa.network_kde(
    network=network,
    points=events,
    bandwidth=500.0,
    lixel_length=50.0,
    method=nsa.NKDVMethod.LION,
    kernel=nsa.NKDVKernel.EPANECHNIKOV,
    normalize=True,
)
~~~

当前公开接口只启用 Epanechnikov 核。ADA 和 LION 是目前插件重点使用的优化
方法，其它可调用方法包括 BASIC、SHARED_SPD、IA 和 ADA_IA。不同实现应在
浮点容差内得到一致趋势，但发布分析前仍建议对目标数据进行回归比较。

lixel_length 越小，输出线段越细，计算量和结果行数也越大。normalize=True
按照每个 lixel 的实际长度计算网络积分并归一化，不是简单除以数组总和。

事件权重示例：

~~~python
kde = nsa.network_kde(
    network=network,
    points=events,
    bandwidth=500.0,
    lixel_length=50.0,
    method="ADA",
    weights=points["severity"].to_numpy(),
)
~~~

权重必须与点数一致、有限且非负，并且至少有一个正值。

## Network K-function

NS/ASPS 与 STARS 方法的原始论文见
[Network K 参考文献](references.md#network-k-function)。

~~~python
k_result = nsa.network_k_function(
    network=network,
    points=events,
    r_min=0.0,
    r_max=5000.0,
    r_step=100.0,
    method=nsa.KFunctionMethod.STARS,
    simulations=99,
    alpha=0.05,
    envelope="global",
    random_state=2026,
)
~~~

可选方法为 SPS、NS 和 STARS。STARS 通常作为较大路网和多阈值任务的优化
选择，但最终仍应以目标数据测试为准。

k_theo 是 Monte Carlo 模拟均值，而不是自动代入的解析理论曲线。因此
simulations=0 时，k_theo、k_lo 和 k_hi 均不存在。

pointwise 包络逐阈值计算分位数，适合查看局部偏离；global 包络基于相对
模拟均值的最大绝对偏差，同时给出全局模拟 p-value，更适合对整条曲线或
整个 surface 作联合判断。

## Spatiotemporal Network K-function

EAR、MTS 和 EARTH 来源于 KDD 2026 EARTH 论文，详见
[ST K 参考文献](references.md#spatiotemporal-network-k-function)。

首先为事件提供数值时间：

~~~python
events = nsa.PointEvents(
    points,
    target_network=network,
    time_col="event_time",
    time_unit="seconds",
    time_origin="Unix epoch (1970-01-01T00:00:00Z)",
    time_mode="datetime",
)
~~~

非数值 datetime 字段会按 UTC 解析并转换为 Unix 秒。数值字段不会自动判断
秒、毫秒、天等单位，调用者必须保证 t 阈值与输入时间使用同一个单位。

~~~python
st_result = nsa.spatiotemporal_network_k_function(
    network=network,
    points=events,
    s_values=[500.0, 1000.0, 1500.0],
    t_values=[86400.0, 172800.0, 259200.0],
    method=nsa.SpatiotemporalKFunctionMethod.EARTH,
    simulations=19,
    random_state=2026,
)
~~~

输出 surface 的 shape 为 (空间阈值数量, 时间阈值数量)。RQS、SPS、EAR、
MTS 和 EARTH 均可选；MTS 与 EARTH 要求 s 和 t 每个轴至少有两个值且等间隔。
时间必须有正跨度，全部相同的时间无法进行 ST K-function。

安装 xarray 后可导出二维结果：

~~~python
surface = st_result.to_xarray("k_obs")
~~~

## 结果保存与恢复

~~~python
result.to_csv("result.csv")
result.to_dataframe()
result.save_npz("result.npz")
restored = nsa.load_result("result.npz")
~~~

CSV 适合交换和人工查看，但不会保存全部元数据与模拟数组；NPZ 用于完整保存
FastNSA 结果类型、数组和统计元数据。不要手动依赖 NPZ 内部键名，应统一使用
load_result 或各结果类的 load_npz。

建议归档道路签名、CRS、输入点来源、映射参数、算法参数、模拟次数、随机种子、
alpha、包络类型、FastNSA 版本，以及 ST 分析的时间单位和原点。
正式发表结果时，还应引用所选算法对应的原始论文，并同时注明 FastNSA 版本。

## 进度与取消

~~~python
cancel_requested = False

def progress(stage, completed, total):
    print(stage, completed, total)
    return cancel_requested

result = nsa.network_k_function(
    network=network,
    points=events,
    progress_callback=progress,
)
~~~

回调返回真值后会请求协作式取消。C++ 后端只会在安全检查点响应，因此长步骤
内部不一定立即停止。取消成功会抛出 FastNSACancelledError。GUI 程序不要在
原生回调里直接修改控件，应把更新转发到主事件循环。

## 常见问题

- 距离明显不合理：检查是否直接使用了经纬度 CRS。
- edge_id 对不上：检查道路是否被排序、过滤、编辑或重新构建。
- 没有 theoretical/envelope：确认 simulations 大于 0。
- ST 提示时间跨度为零：检查时间解析、单位和字段内容。
- 导入原生库失败：运行 python -m fastnsa._selftest，并确认 wheel 与 Python
  的系统、CPU 架构和位数一致。
- 运行很慢：复用道路和 mapping，减少探索阶段的阈值网格与模拟次数，并使用
  Release 构建。

更多排错步骤见 [Troubleshooting](troubleshooting.md)。

