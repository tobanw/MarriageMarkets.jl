# verify that super/sub-modularity of production implies positive/negative assortative matching.

# get results from worker processes
pam = fetch(pam_job)
nam = fetch(nam_job)

# sum of diag corners greater than sum of anti-diag corners
@fact pam.matches[1,1] + pam.matches[end,end] -->
      greater_than(pam.matches[1,end] + pam.matches[end,1]) "expected positive assortativity"
@fact nam.matches[1,1] + nam.matches[end,end] -->
      less_than(nam.matches[1,end] + nam.matches[end,1]) "expected negative assortativity"

# symmetry
@fact pam.matches --> roughly(pam.matches') "expected symmetry of matches"
@fact pam.msingle --> roughly(pam.fsingle) "expected symmetry of singles"
