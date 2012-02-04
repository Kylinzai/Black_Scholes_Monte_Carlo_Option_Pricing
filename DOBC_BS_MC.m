% Ioannis Anagnostou, 2010

% function DOBC_BS_MC uses Monte Carlo simulation to estimate 
% the price of a Down-and-out Barrier Call Option under the
% Black Scholes Model

% S0: initial price of the underlying
% K: option strike
% H: barrier
% r: annualized risk-free interest rate
% q: divident
% sigma: volatility of the underlying
% T: time to maturity
% N: time steps
% n: number of simulations

function DOBC_BS_MC(S0,K,H,r,q,sigma,T,N,n)

dt=T/N;
Z=randn(n,N);
dW=sqrt(dt)*Z;
S=S0*cumprod(exp((r-q-0.5*(sigma^2))*dt+sigma*dW),2);
ST=S(:,N);
m=min(S,[],2); 

for i=1:n  
    if m(i)>=H
        c(i)=max(0,ST(i)-K);
    else 
        c(i)=0;
    end
end 
    
call_price=exp(-r*T)*mean(c)