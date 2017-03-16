# uses variables from setup_tests.jl

srmgmkt = fetch(srsymm_multi_job)

# check symmetry
@fact srmgmkt.u_m --> roughly(srmgmkt.u_f) "expected symmetry of singles"

# check that supermodular surplus results in positive assortative matching
@fact srmgmkt.α[1,1,1,1] + srmgmkt.α[end,end,end,end] -->
      greater_than(srmgmkt.α[1,1,end,end] + srmgmkt.α[end,end,1,1]) "expected positive assortativity"

