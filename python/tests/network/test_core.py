def test_network_ir_structure(simple_grid_network):
    net = simple_grid_network
    # 检查节点重映射是否正确
    assert net.num_nodes == 4
    # 检查 CSR 偏移量：4个节点，offsets 长度应为 5
    assert len(net._data_cache['offsets']) == 5
    # 检查总长度是否正确 (100 + 100)
    assert net.total_length == 200.0
    