library(UTILS)
library(ZOOPS)
k = 14
fasta = UTILS::generate_ZOOPS(
  clevel = 0.95, 
  k = k, 
  n = 5000, 
  t = 150,
  pos = 50,
  noise = 0.3
)

n = length(fasta)
t = nchar(fasta[1])
m = t - k + 1

alpha = t(extraDistr::rdirichlet(k, c(1,1,1,1)))
beta = createMarkovChain(fasta)

new_alpha = zoops(fasta, alpha, beta, 1e-2, 1e2, 0.5)
UTILS::plot_data(new_alpha)

