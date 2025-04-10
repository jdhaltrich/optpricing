# option pricing functions

black_scholes <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
    black_scholes_call <- S*pnorm(d1,0,1,lower.tail = TRUE)-K*pnorm(d2,0,1,lower.tail = TRUE)*exp(-r*delta_t)
    black_scholes_put <- K*pnorm(-d2,0,1,lower.tail = TRUE)*exp(-r*delta_t)-S*pnorm(-d1,0,1,lower.tail = TRUE)
  list(black_scholes_call,black_scholes_put)
}

greek_delta <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  delta_call <- pnorm(d1,0,1,lower.tail = TRUE)
  delta_put <- -1*pnorm(-d1,0,1,lower.tail = TRUE)
  list(delta_call,delta_put)
}

greek_gamma <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
  (1/(S*sigma*(delta_t^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))
}

greek_theta <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
  theta_call <- -S*(1/(2*(delta_t^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))*sigma-r*K*exp(-r*delta_t)*pnorm(d2,0,1,lower.tail = TRUE)
  theta_put <- -S*(1/(2*(delta_t^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))*sigma+r*K*exp(-r*delta_t)*pnorm(-d2,0,1,lower.tail = TRUE)
  list(theta_call,theta_put)
}

greek_vega <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
    S*(delta_t^(.5))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))
}

greek_rho <- function(S,sigma,delta_t,K,r){
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
    rho_call <- K*delta_t*exp(-.5*delta_t)*pnorm(d2,0,1,lower.tail = TRUE)
    rho_put <- -1*K*delta_t*exp(-.5*delta_t)*pnorm(-d2,0,1,lower.tail = TRUE)
    list(rho_call,rho_put)
}

implied_volatility <- function(init_vol,delta_t,S,K,r,call_real) {
  iteration_1 <- init_vol
  d1_1 <- (1/(iteration_1*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_1^2))*delta_t)
  d2_1 <- (1/(iteration_1*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_1^2))*delta_t)
  black_scholes_call_iteration_1 <- S*pnorm(d1_1,0,1,lower.tail = TRUE)-K*pnorm(d2_1,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_1 <- black_scholes_call_iteration_1-call_real
  f_1_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_1^2)*(delta_t^(0.5))
  newton_raphson_1 <- iteration_1 - f_1/f_1_derivative
  
  iteration_2 <- newton_raphson_1
  d1_2 <- (1/(iteration_2*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_2^2))*delta_t)
  d2_2 <- (1/(iteration_2*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_2^2))*delta_t)
  black_scholes_call_iteration_2 <- S*pnorm(d1_2,0,1,lower.tail = TRUE)-K*pnorm(d2_2,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_2 <- black_scholes_call_iteration_2-call_real
  f_2_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_2^2)*(delta_t^(0.5))
  newton_raphson_2 <- iteration_2 - f_2/f_2_derivative
  
  iteration_3 <- newton_raphson_2
  d1_3 <- (1/(iteration_3*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_3^2))*delta_t)
  d2_3 <- (1/(iteration_3*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_3^2))*delta_t)
  black_scholes_call_iteration_3 <- S*pnorm(d1_3,0,1,lower.tail = TRUE)-K*pnorm(d2_3,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_3 <- black_scholes_call_iteration_3-call_real
  f_3_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_3^2)*(delta_t^(0.5))
  newton_raphson_3 <- iteration_3 - f_3/f_3_derivative
  
  iteration_4 <- newton_raphson_3
  d1_4 <- (1/(iteration_4*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_4^2))*delta_t)
  d2_4 <- (1/(iteration_4*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_4^2))*delta_t)
  black_scholes_call_iteration_4 <- S*pnorm(d1_4,0,1,lower.tail = TRUE)-K*pnorm(d2_4,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_4 <- black_scholes_call_iteration_4-call_real
  f_4_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_4^2)*(delta_t^(0.5))
  newton_raphson_4 <- iteration_4 - f_4/f_4_derivative

  iteration_5 <- newton_raphson_4
  d1_5 <- (1/(iteration_5*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_5^2))*delta_t)
  d2_5 <- (1/(iteration_5*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_5^2))*delta_t)
  black_scholes_call_iteration_5 <- S*pnorm(d1_5,0,1,lower.tail = TRUE)-K*pnorm(d2_5,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_5 <- black_scholes_call_iteration_5-call_real
  f_5_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_5^2)*(delta_t^(0.5))
  newton_raphson_5 <- iteration_5 - f_5/f_5_derivative
  
  iteration_6 <- newton_raphson_5
  d1_6 <- (1/(iteration_6*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_6^2))*delta_t)
  d2_6 <- (1/(iteration_6*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_6^2))*delta_t)
  black_scholes_call_iteration_6 <- S*pnorm(d1_6,0,1,lower.tail = TRUE)-K*pnorm(d2_6,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_6 <- black_scholes_call_iteration_6-call_real
  f_6_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_6^2)*(delta_t^(0.5))
  newton_raphson_6 <- iteration_6 - f_6/f_6_derivative

  iteration_7 <- newton_raphson_6
  d1_7 <- (1/(iteration_7*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_7^2))*delta_t)
  d2_7 <- (1/(iteration_7*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_7^2))*delta_t)
  black_scholes_call_iteration_7 <- S*pnorm(d1_7,0,1,lower.tail = TRUE)-K*pnorm(d2_7,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_7 <- black_scholes_call_iteration_7-call_real
  f_7_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_7^2)*(delta_t^(0.5))
  newton_raphson_7 <- iteration_7 - f_7/f_7_derivative

  iteration_8 <- newton_raphson_7
  d1_8 <- (1/(iteration_8*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_8^2))*delta_t)
  d2_8 <- (1/(iteration_8*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_8^2))*delta_t)
  black_scholes_call_iteration_8 <- S*pnorm(d1_8,0,1,lower.tail = TRUE)-K*pnorm(d2_8,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_8 <- black_scholes_call_iteration_8-call_real
  f_8_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_8^2)*(delta_t^(0.5))
  newton_raphson_8 <- iteration_8 - f_8/f_8_derivative

  iteration_9 <- newton_raphson_8
  d1_9 <- (1/(iteration_9*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_9^2))*delta_t)
  d2_9 <- (1/(iteration_9*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_9^2))*delta_t)
  black_scholes_call_iteration_9 <- S*pnorm(d1_9,0,1,lower.tail = TRUE)-K*pnorm(d2_9,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_9 <- black_scholes_call_iteration_9-call_real
  f_9_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_9^2)*(delta_t^(0.5))
  newton_raphson_9 <- iteration_9 - f_9/f_9_derivative

  iteration_10 <- newton_raphson_9
  d1_10 <- (1/(iteration_10*(delta_t)^(.5)))*(log(S/K)+(r+.5*(iteration_10^2))*delta_t)
  d2_10 <- (1/(iteration_10*(delta_t)^(.5)))*(log(S/K)+(r-.5*(iteration_10^2))*delta_t)
  black_scholes_call_iteration_10 <- S*pnorm(d1_10,0,1,lower.tail = TRUE)-K*pnorm(d2_10,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_10 <- black_scholes_call_iteration_10-call_real
  f_10_derivative <- S*(1/((2*pi)^(.5)))*exp(-.5*d1_10^2)*(delta_t^(0.5))
  newton_raphson_10 <- iteration_10 - f_10/f_10_derivative  
  
  c(newton_raphson_1,newton_raphson_2,newton_raphson_3,newton_raphson_4,newton_raphson_5,newton_raphson_6,newton_raphson_7,newton_raphson_8,
    newton_raphson_9,newton_raphson_10)
  
}
