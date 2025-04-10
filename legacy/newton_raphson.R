library(RMySQL)

implied_volatility <- function(S,sigma,delta_t,K,r,call_real) {
  iteration_1 <- S 
  d1_1 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_1/K)+(r+.5*(sigma^2))*delta_t)
  d2_1 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_1/K)+(r-.5*(sigma^2))*delta_t)
  black_scholes_call_it_1 <- iteration_1*pnorm(d1_1,0,1,lower.tail = TRUE)-K*pnorm(d2_1,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_1 <- black_scholes_call_it_1-call_real
  f_1_derivative <- iteration_1*(1/((2*pi)^(.5)))*exp(-.5*d1_1^2)*(delta_t^(0.5))
  newton_raphson_1 <- iteration_1 - f_1/f_1_derivative
 
  iteration_2 <- newton_raphson_1 
  d1_2 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_2/K)+(r+.5*(sigma^2))*delta_t)
  d2_2 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_2/K)+(r-.5*(sigma^2))*delta_t)
  black_scholes_call_it_2 <- iteration_2*pnorm(d1_2,0,1,lower.tail = TRUE)-K*pnorm(d2_2,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_2 <- black_scholes_call_it_2-call_real
  f_2_derivative <- iteration_2*(1/((2*pi)^(.5)))*exp(-.5*d1_2^2)*(delta_t^(0.5))
  newton_raphson_2 <- iteration_2 - f_2/f_2_derivative
 
  iteration_3 <- newton_raphson_2 
  d1_3 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_3/K)+(r+.5*(sigma^2))*delta_t)
  d2_3 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_3/K)+(r-.5*(sigma^2))*delta_t)
  black_scholes_call_it_3 <- iteration_3*pnorm(d1_3,0,1,lower.tail = TRUE)-K*pnorm(d2_3,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_3 <- black_scholes_call_it_3-call_real
  f_3_derivative <- iteration_3*(1/((2*pi)^(.5)))*exp(-.5*d1_3^2)*(delta_t^(0.5))
  newton_raphson_3 <- iteration_3 - f_3/f_3_derivative
 
  iteration_4 <- newton_raphson_3 
  d1_4 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_4/K)+(r+.5*(sigma^2))*delta_t)
  d2_4 <- (1/(sigma*(delta_t)^(.5)))*(log(iteration_4/K)+(r-.5*(sigma^2))*delta_t)
  black_scholes_call_it_4 <- iteration_4*pnorm(d1_4,0,1,lower.tail = TRUE)-K*pnorm(d2_4,0,1,lower.tail = TRUE)*exp(-r*delta_t)
  f_4 <- black_scholes_call_it_4-call_real
  f_4_derivative <- iteration_4*(1/((2*pi)^(.5)))*exp(-.5*d1_4^2)*(delta_t^(0.5))
  newton_raphson_4 <- iteration_4 - f_4/f_4_derivative

  
  c(newton_raphson_1,newton_raphson_2,newton_raphson_3,newton_raphson_4)
}

