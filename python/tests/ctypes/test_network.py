def test_network_create_destroy(simple_network):
    # 只要 fixture 没 segfault，这个 test 就算通过
    assert simple_network is not None
