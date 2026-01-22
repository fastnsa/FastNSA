#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <cmath>
#include <iomanip>

#include "fastnsa_capi.h"

//  --HELPER FUNCTION--
// 格式: [ByteOrder:1][Type:4][NumPoints:4][X:8][Y:8]...



std::vector<uint8_t> create_linestring_wkb(double x1, double y1, double x2, double y2) {
    std::vector<uint8_t> wkb;
    
    // 1. Byte Order (1 = Little Endian)
    wkb.push_back(1);
    
    // 2. Type (2 = LineString)
    uint32_t type = 2;
    uint8_t* p_type = reinterpret_cast<uint8_t*>(&type);
    wkb.insert(wkb.end(), p_type, p_type + 4);
    
    // 3. Num Points (2 points)
    uint32_t num_points = 2;
    uint8_t* p_num = reinterpret_cast<uint8_t*>(&num_points);
    wkb.insert(wkb.end(), p_num, p_num + 4);
    
    // 4. Points (Double x, Double y)
    double coords[] = {x1, y1, x2, y2};
    uint8_t* p_coords = reinterpret_cast<uint8_t*>(coords);
    wkb.insert(wkb.end(), p_coords, p_coords + sizeof(coords));
    
    return wkb;
}


int main() {
    std::cout << "=== 开始 FastNSA C-API 测试 ===" << std::endl;

    // -------------------------------------------------
    // 1. 准备拓扑数据
    // 网络结构: Node 0 --(edge0)--> Node 1 --(edge1)--> Node 2
    // 坐标示意: (0,0) ----------- (10,0) ----------- (20,0)
    // -------------------------------------------------
    int n_nodes = 3;
    int n_edges = 2;

    // Edges input format: [u, v, length]
    std::vector<double> edges = {
        0, 1, 10.0,  // Edge 0
        1, 2, 10.0   // Edge 1
    };

    // CSR (Adjacency)
    // Node 0 -> edges to 1 (let's assume index 0)
    // Node 1 -> edges to 2 (let's assume index 1)
    // Node 2 -> no outgoing
    std::vector<int> csr_offsets = {0, 1, 2, 2}; // size = n_nodes + 1
    std::vector<int> csr_indices = {1, 2};       // Targets

    std::cout << "[Step 1] 创建网络..." << std::endl;
    void* net = fastnsa_network_create(edges.data(), n_edges, csr_offsets.data(), csr_indices.data(), n_nodes);
    assert(net != nullptr);

    // -------------------------------------------------
    // 2. 设置几何数据 (WKB)
    // 这一步对 fastnsa_map_points 至关重要
    // -------------------------------------------------
    std::cout << "[Step 2] 设置 WKB..." << std::endl;
    
    // 生成 Edge 0 的几何: (0,0) -> (10,0)
    auto wkb0 = create_linestring_wkb(0.0, 0.0, 10.0, 0.0);
    // 生成 Edge 1 的几何: (10,0) -> (20,0)
    auto wkb1 = create_linestring_wkb(10.0, 0.0, 20.0, 0.0);

    // 拼接所有 WKB 数据
    std::vector<uint8_t> all_wkb_data;
    all_wkb_data.insert(all_wkb_data.end(), wkb0.begin(), wkb0.end());
    all_wkb_data.insert(all_wkb_data.end(), wkb1.begin(), wkb1.end());

    // 计算 offsets (byte positions)
    std::vector<int> wkb_offsets = {
        0, 
        (int)wkb0.size(), 
        (int)(wkb0.size() + wkb1.size())
    };

    fastnsa_network_set_wkbs(net, all_wkb_data.data(), wkb_offsets.data(), n_edges); //TODO 

    // -------------------------------------------------
    // 3. 测试点映射 (Map Points)
    // 点 P1: (5, 1) -> 应该投影到 Edge 0, offset 5.0
    // 点 P2: (15, -1) -> 应该投影到 Edge 1, offset 5.0
    // -------------------------------------------------
    std::cout << "[Step 3] 映射点..." << std::endl;
    std::vector<double> px = {5.0, 15.0};
    std::vector<double> py = {1.0, -1.0};
    int n_points = 2;

    void* pts = fastnsa_map_points(net, px.data(), py.data(), n_points);
    std::cout << "[Step 3] 映射点完成" << std::endl;  //test
    assert(pts != nullptr);
    
    // 注意: 这里无法直接验证 pts 内部数据，因为是指针。
    // 如果你有 getter 函数可以验证 edge_id 和 offset 是否正确。
    // 或者通过后续计算 K 函数没有崩溃来间接验证。

    // -------------------------------------------------
    // 4. (可选) 测试手动设置点集
    // -------------------------------------------------
    std::cout << "[Step 4] 手动设置点集测试..." << std::endl;
    void* manual_pts = fastnsa_pointset_create(net);
    std::vector<int> m_ids = {0, 1};
    std::vector<double> m_offsets = {5.0, 5.0};
    
    int set_res = fastnsa_pointset_set(manual_pts, m_ids.data(), m_offsets.data(), 2);
    if(set_res == 0) std::cout << "  -> 设置成功" << std::endl;
    else std::cout << "  -> 设置失败 Code: " << set_res << std::endl;
    assert(set_res == 0);

    // 测试越界
    std::vector<double> bad_offsets = {15.0, 5.0}; // Edge 0 只有 10.0 长
    assert(fastnsa_pointset_set(manual_pts, m_ids.data(), bad_offsets.data(), 2) == -2);
    std::cout << "  -> 越界检查通过" << std::endl;

    // -------------------------------------------------
    // 5. 计算 K 函数
    // -------------------------------------------------
    std::cout << "[Step 5] 计算 K 函数..." << std::endl;
    std::vector<double> r_values = {1.0, 2.0, 5.0, 10.0};
    int n_r = r_values.size();
    std::vector<double> out_k(n_r);
    
    // 使用自动映射的点集 pts 进行计算
    fastnsa_compute_k_function(
        net, 
        pts, 
        r_values.data(), 
        n_r, 
        FASTNSA_K_NS, 
        0, // No simulations
        out_k.data(), 
        nullptr
    );

    std::cout << "  -> K Values: ";
    for (double v : out_k) std::cout << std::fixed << std::setprecision(4) << v << " ";
    std::cout << std::endl;

    // -------------------------------------------------
    // 6. 销毁内存
    // -------------------------------------------------
    std::cout << "[Step 6] 清理..." << std::endl;
    fastnsa_pointset_destroy(pts);
    fastnsa_pointset_destroy(manual_pts);
    fastnsa_network_destroy(net);

    std::cout << "=== 测试完成 ===" << std::endl;
    return 0;
}

