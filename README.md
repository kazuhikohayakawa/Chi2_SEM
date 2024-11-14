# chi2test_SEM

test for $\chi^2$
In a structural equation modeling, a model-implied covariance matrix $\boldmath{\Sigma}(\boldmath{\theta}), (p \times p)$ is fitted to the sample covariance matrix $\bm{S}_N$ which is obtained from $N$ observations $(\bm{z}_1,...,\bm{z}_N),(p \times 1)$ in order to minimize a discrepancy function. In practice, unknown parameter $\bmg{\theta}$ is typically obtained by minimizing the Wishart likelihood $F_{ML}(\bmg{\theta})$: 
\begin{eqnarray*}
\widehat{\bmg{\theta}}_{ML} &=& \mathop{\argmin}_{\theta }F_{ML}\left( \bmg{\theta} \right), \label{eq_ML} \\
F_{ML} \left( \bmg{\theta} \right) &=& \log \left\vert \bmg{\Sigma}\left( \bmg{\theta}\right) \right\vert - \bm{\log }\left\vert \bm{S}_N \right\vert + \tr \left( \bm{S}_N\bmg{\Sigma}\left( \bmg{\theta}\right)^{-1}\right) -p. \label{eq_F_ML}
\end{eqnarray*}
The conventional goodness of fit test or a likelihood ratio test is defined as 
\begin{eqnarray}
T_{ML} = n \cdot \widehat{F}_{ML}  
\end{eqnarray}
where $\widehat{F}_{ML}=F_{ML}( \widehat{\bmg{\theta}}_{ML} )$ and $n=N-1$.
