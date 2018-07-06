
__global__ void mykernel( double *Gamma, double *g, double *Theta, double *X, int N, int K, int T) {

	/* compute kernel index */
	int t = blockIdx.x*blockDim.x + threadIdx.x;

	if(t<T){
		int mink;
		
		/* compute g(:,t) */
		for(int k=0;k<K;k++){
			/* compute dot product g(k,t) = <X(:,t) - Theta(:,k),X(:,t) - Theta(:,k)> */
			g[t*K+k] = 0;
			for(int n=0;n<N;n++){
				g[t*K+k] += (X[t*N+n] - Theta[k*N+n])*(X[t*N+n] - Theta[k*N+n]);
			}
			
			/* if this is first row, then Gamma(k,t) is minimal value */
			if(k==0){
				mink=0; /* index k with min value of g(:,t) */
				Gamma[t*K+k] = 1;
			} else {
				/* is this smaller value then previous one? */
				if(g[t*K+k] < g[t*K+mink]){
					/* old one is not min, set it equal to zero */
					Gamma[t*K+mink] = 0;
					mink=k;
					Gamma[t*K+k] = 1;
				} else {
					/* it is not min */
					Gamma[t*K+k] = 0;
				}
			}
		}

	}
}
