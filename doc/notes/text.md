---
title: Notes
author: Sinan Shi 
date: \today 
csl: bioinformatics.csl
---

# EM Algorithm

## The General EM Algorithm
Given a joint distribution $p(\mathbf{X,Z|\theta})$ over observed variables $\mathbf{X}$ and latent variables $\mathbf{Z}$, 
governed by parameters $\mathbf{\theta}$, the goal is to maximize the likelihood function $p\mathbf{(X|\theta)}$.

1. Choose an initial setting for the parameters $\theta^{old}$.
2. **E step** Evaluate $p(\mathbf{Z|X, \theta^{old}})$
3. **M step** Evaluate $\mathbf{\theta^{new}}$ given by
$$\mathbf{\theta^{new}} = \argmax_x Q(\mathbf{\theta, \theta^{old}})$$

$$ Q(\mathbf{\theta}, \mathbf{\theta}^{old}) = 
\sum_{\mathbf{z}} p(\mathbf{Z}|\mathbf{X}, \mathbf{\theta}^{old}) \log p(\mathbf{X},\mathbf{Z}|\mathbf{\theta})$$ {#eq:q-function}


## HMM

$$p(\mathbf{X,Z}|\theta) = p(\mathbf{z}_1|\mathbf{\pi}) \prod_{n=2}^{N}p(\mathbf{z}_n|\mathbf{z}_{n-1}) \prod_{m=1}^{N}p(\mathbf{x}_m|\mathbf{z}_{m}, \mathbf{\phi})$$ {#eq:hmm-complete-likelihood}

$$\gamma(\mathbf{z}_n) = p(\mathbf{z}_n|\mathbf{X}, \mathbf{\theta}^{old})$$

$$\xi(\mathbf{z}_{n-1}, \mathbf{z}_n) = p(\mathbf{z}_{n-1}, \mathbf{z}_n|\mathbf{X}, \theta^{old})$$

[@em_tutorial]
[@Bishop]

# Reference
