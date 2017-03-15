# uses variables from setup_tests.jl

mgmkt = fetch(symm_multi_job)

# check symmetry
@fact mgmkt.msingle --> roughly(mgmkt.fsingle) "expected symmetry of singles"

mgmkt2 = fetch(asym_multi_job)

# check that supermodular surplus results in positive assortative matching
@fact mgmkt2.matches[1,1,1,1] + mgmkt2.matches[end,end,end,end] -->
      greater_than(mgmkt2.matches[1,1,end,end] + mgmkt2.matches[end,end,1,1]) "expected positive assortativity"
