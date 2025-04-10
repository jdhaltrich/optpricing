library(RMySQL)
library(xlsx)

option_database <- dbConnect(MySQL(),user='jdh',password='brha1j9d',dbname='options',host='127.0.0.1')

instrument <- dbGetQuery(option_database,
                         "call options.instrument_analysis(20180101,20200330,'AAPL200320C00500000','nasdaq',20180101,20200330,'AAPL.OQ')"
)

instrument_2 <- dbGetQuery(option_database,
                           "
                           select *
                           from options.extract_instrument_analysis
                           where effective_date<=20200211 and ric is not null
                           
                           union all
                           
                           select *
                           from options.extract_instrument_analysis
                           where effective_date>20200211 and flag_no_price_cboe=1
                           
                           order by effective_date
                           
                           "
)

stock_price <- instrument_2[,"precio_cierre"]

prices_length <- length(stock_price)

parameter_days_mean <- 20

stock_prices <- lapply(
  (parameter_days_mean+1):prices_length,
  function(i) {stock_price[(i-parameter_days_mean):i]}
)

stock_prices_2 <- matrix(
  unlist(stock_prices),
  ncol=prices_length-21+1,
  byrow=FALSE
)

stock_prices_lag <- stock_prices_2[1:(nrow(stock_prices_2)-1),,drop=FALSE]

stock_prices_2 <- stock_prices_2[-1,,drop=FALSE]

s <- tail(
  ifelse(is.na(stock_price),0,stock_price),
  prices_length-parameter_days_mean
)

log_variation <- log(stock_prices_2*(1/stock_prices_lag))

mean_log_variation_suma <- apply(
  log_variation,
  2,
  sum
)

mean_log_variation_length <- apply(
  log_variation,
  2,
  length
)

mean_log_variation <- mean_log_variation_suma/mean_log_variation_length

mean_matrix <- matrix(
  unlist(lapply(
    1:(prices_length-parameter_days_mean),
    function(i) {rep(mean_log_variation[i],20)})
  ),
  ncol=(prices_length-parameter_days_mean),
  byrow = FALSE
)

mean_substraction <- log_variation-mean_matrix

mean_sq <- mean_substraction^2

sum_sq_1 <- apply(
  mean_sq,
  2,
  sum
)

sum_sq_2 <- apply(
  mean_sq,
  2,
  length
)

sum_sq_3 <- sum_sq_1/sum_sq_2

volatility <- ((sum_sq_3)^(.5))*(252^(.5))

flag_remaining_days <- ifelse(
  is.na(instrument_2[,"RIC"]),
  1,
  0
)

remaining_days <- sum(flag_remaining_days)

strike_final <- max(
  ifelse(
    is.na(instrument_2[,"strike"]),
    0,
    instrument_2[,"strike"]
  )
)

expiration_date <- max(
  ifelse(
    is.na(instrument_2[,"expiration_date"]),
    0,
    instrument_2[,"expiration_date"]*1
  )
)

flag_live_option <- tail(
  ifelse(
    (instrument_2[,"effective_date"])<=expiration_date,
    1,
    0
  ),
  prices_length-parameter_days_mean
)

delta_t_year <- unlist(
  lapply(
    1:(prices_length-parameter_days_mean),
    function(i) {sum(flag_live_option[i:length(flag_live_option)])})
)*(1/252)

black_scholes <- function(s,volatility,delta_t_year,k,r) {
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  d2 <- d1 - volatility*(delta_t_year^(.5))
  black_scholes_call<-s*pnorm(d1,0,1,lower.tail = TRUE)-k*pnorm(d2,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  black_scholes_put<-k*pnorm(-d2,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)-s*pnorm(-d1,0,1,lower.tail = TRUE)
  list(black_scholes_call,black_scholes_put)
}

greek_delta <- function(s,volatility,delta_t_year,k,r) {
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  delta_call <- pnorm(d1,0,1,lower.tail = TRUE)
  delta_put <- -1*pnorm(-d1,0,1,lower.tail = TRUE)
  list(delta_call,delta_put)
}

greek_gamma <- function(s,volatility,delta_t_year,k,r) {
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  d2 <- d1 - volatility*(delta_t_year^(.5))
  (1/(s*volatility*(delta_t_year^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))
}

greek_theta <- function(s,volatility,delta_t_year,k,r) {
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  d2 <- d1 - volatility*(delta_t_year^(.5))
  theta_call <- -s*(1/(2*(delta_t_year^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))*volatility-r*k*exp(-r*delta_t_year)*pnorm(d2,0,1,lower.tail = TRUE)
  theta_put <- -s*(1/(2*(delta_t_year^(.5))))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))*volatility+r*k*exp(-r*delta_t_year)*pnorm(-d2,0,1,lower.tail = TRUE)
  list(theta_call,theta_put)
}

greek_vega <- function(s,volatility,delta_t_year,k,r) {
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  d2 <- d1 - volatility*(delta_t_year^(.5))
  s*(delta_t_year^(.5))*(1/((2*pi)^(.5)))*exp(-.5*(d1^2))
}

greek_rho <- function(s,volatility,delta_t_year,k,r){
  d1 <- (1/(volatility*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(volatility^2))*delta_t_year)
  d2 <- d1 - volatility*(delta_t_year^(.5))
  rho_call <- k*delta_t_year*exp(-.5*delta_t_year)*pnorm(d2,0,1,lower.tail = TRUE)
  rho_put <- -1*k*delta_t_year*exp(-.5*delta_t_year)*pnorm(-d2,0,1,lower.tail = TRUE)
  list(rho_call,rho_put)
}

black_scholes_strip <- mapply(
  black_scholes,
  s,volatility,delta_t_year,strike_final,.03
)

black_scholes_prices_matrix <- matrix(
  unlist(black_scholes_strip),
  ncol=(prices_length-parameter_days_mean)
)

calls_prices <- black_scholes_prices_matrix[1,]

puts_prices <- black_scholes_prices_matrix[2,]


greek_delta_strip <- mapply(
  greek_delta,
  s,volatility,delta_t_year,strike_final,.03
)

greek_delta_matrix <- matrix(
  unlist(greek_delta_strip),
  ncol=(prices_length-parameter_days_mean)
)

greek_delta_call <- greek_delta_matrix[1,]

greek_delta_put <- greek_delta_matrix[2,]


greek_gamma_strip <- mapply(
  greek_gamma,
  s,volatility,delta_t_year,strike_final,.03
)

greek_theta_strip <- mapply(
  greek_theta,
  s,volatility,delta_t_year,strike_final,.03
)

greek_theta_matrix <- matrix(
  unlist(greek_theta_strip),
  ncol=(prices_length-parameter_days_mean)
)

greek_theta_call <- greek_theta_matrix[1,]

greek_theta_put <- greek_theta_matrix[2,]


greek_vega_strip <- mapply(
  greek_vega,
  s,volatility,delta_t_year,strike_final,.03
)


greek_rho_strip <- mapply(
  greek_rho,
  s,volatility,delta_t_year,strike_final,.03
)

greek_rho_matrix <- matrix(
  unlist(greek_rho_strip),
  ncol=(prices_length-parameter_days_mean)
)

greek_rho_call <- greek_rho_matrix[1,]

greek_rho_put <- greek_rho_matrix[2,]


implied_volatility <- function(input_volatility,delta_t_year,s,k,r,call_real) {
  initial_guess <- input_volatility
  d1_1 <- (1/(initial_guess*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(initial_guess^2))*delta_t_year)
  d2_1 <- (1/(initial_guess*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(initial_guess^2))*delta_t_year)
  black_scholes_call_init_guess <- s*pnorm(d1_1,0,1,lower.tail = TRUE)-k*pnorm(d2_1,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_1 <- black_scholes_call_init_guess-call_real
  f_1_derivative_1 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_1^2)*(delta_t_year^(0.5))
  newton_raphson_1 <- initial_guess - f_1/f_1_derivative_1
  
  guess_2 <- newton_raphson_1
  d1_2 <- (1/(guess_2*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_2^2))*delta_t_year)
  d2_2 <- (1/(guess_2*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_2^2))*delta_t_year)
  black_scholes_call_guess_2 <- s*pnorm(d1_2,0,1,lower.tail = TRUE)-k*pnorm(d2_2,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_2 <- black_scholes_call_guess_2-call_real
  f_1_derivative_2 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_2^2)*(delta_t_year^(0.5))
  newton_raphson_2 <- guess_2 - f_2/f_1_derivative_2
  
  guess_3 <- newton_raphson_2
  d1_3 <- (1/(guess_3*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_3^2))*delta_t_year)
  d2_3 <- (1/(guess_3*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_3^2))*delta_t_year)
  black_scholes_call_guess_3 <- s*pnorm(d1_3,0,1,lower.tail = TRUE)-k*pnorm(d2_3,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_3 <- black_scholes_call_guess_3-call_real
  f_1_derivative_3 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_3^2)*(delta_t_year^(0.5))
  newton_raphson_3 <- guess_3 - f_3/f_1_derivative_3
  
  guess_4 <- newton_raphson_3
  d1_4 <- (1/(guess_4*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_4^2))*delta_t_year)
  d2_4 <- (1/(guess_4*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_4^2))*delta_t_year)
  black_scholes_call_guess_4 <- s*pnorm(d1_4,0,1,lower.tail = TRUE)-k*pnorm(d2_4,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_4 <- black_scholes_call_guess_4-call_real
  f_1_derivative_4 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_4^2)*(delta_t_year^(0.5))
  newton_raphson_4 <- guess_4 - f_4/f_1_derivative_4
  
  guess_5 <- newton_raphson_4
  d1_5 <- (1/(guess_5*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_5^2))*delta_t_year)
  d2_5 <- (1/(guess_5*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_5^2))*delta_t_year)
  black_scholes_call_guess_5 <- s*pnorm(d1_5,0,1,lower.tail = TRUE)-k*pnorm(d2_5,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_5 <- black_scholes_call_guess_5-call_real
  f_1_derivative_5 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_5^2)*(delta_t_year^(0.5))
  newton_raphson_5 <- guess_5 - f_5/f_1_derivative_5
  
  guess_6 <- newton_raphson_5
  d1_6 <- (1/(guess_6*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_6^2))*delta_t_year)
  d2_6 <- (1/(guess_6*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_6^2))*delta_t_year)
  black_scholes_call_guess_6 <- s*pnorm(d1_6,0,1,lower.tail = TRUE)-k*pnorm(d2_6,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_6 <- black_scholes_call_guess_6-call_real
  f_1_derivative_6 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_6^2)*(delta_t_year^(0.5))
  newton_raphson_6 <- guess_6 - f_6/f_1_derivative_6
  
  guess_7 <- newton_raphson_6
  d1_7 <- (1/(guess_7*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_7^2))*delta_t_year)
  d2_7 <- (1/(guess_7*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_7^2))*delta_t_year)
  black_scholes_call_guess_7 <- s*pnorm(d1_7,0,1,lower.tail = TRUE)-k*pnorm(d2_7,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_7 <- black_scholes_call_guess_7-call_real
  f_1_derivative_7 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_7^2)*(delta_t_year^(0.5))
  newton_raphson_7 <- guess_7 - f_7/f_1_derivative_7
  
  guess_8 <- newton_raphson_7
  d1_8 <- (1/(guess_8*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_8^2))*delta_t_year)
  d2_8 <- (1/(guess_8*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_8^2))*delta_t_year)
  black_scholes_call_guess_8 <- s*pnorm(d1_8,0,1,lower.tail = TRUE)-k*pnorm(d2_8,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_8 <- black_scholes_call_guess_8-call_real
  f_1_derivative_8 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_8^2)*(delta_t_year^(0.5))
  newton_raphson_8 <- guess_8 - f_8/f_1_derivative_8
  
  guess_9 <- newton_raphson_8
  d1_9 <- (1/(guess_9*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_9^2))*delta_t_year)
  d2_9 <- (1/(guess_9*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_9^2))*delta_t_year)
  black_scholes_call_guess_9 <- s*pnorm(d1_9,0,1,lower.tail = TRUE)-k*pnorm(d2_9,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_9 <- black_scholes_call_guess_9-call_real
  f_1_derivative_9 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_9^2)*(delta_t_year^(0.5))
  newton_raphson_9 <- guess_9 - f_9/f_1_derivative_9
  
  guess_10 <- newton_raphson_9
  d1_10 <- (1/(guess_10*(delta_t_year)^(.5)))*(log(s/k)+(r+.5*(guess_10^2))*delta_t_year)
  d2_10 <- (1/(guess_10*(delta_t_year)^(.5)))*(log(s/k)+(r-.5*(guess_10^2))*delta_t_year)
  black_scholes_call_guess_10 <- s*pnorm(d1_10,0,1,lower.tail = TRUE)-k*pnorm(d2_10,0,1,lower.tail = TRUE)*exp(-r*delta_t_year)
  f_10 <- black_scholes_call_guess_10-call_real
  f_1_derivative_10 <- s*(1/((2*pi)^(.5)))*exp(-.5*d1_10^2)*(delta_t_year^(0.5))
  newton_raphson_10 <- guess_10 - f_10/f_1_derivative_10
  
  c(newton_raphson_1,newton_raphson_2,newton_raphson_3,newton_raphson_4,newton_raphson_5,newton_raphson_6,newton_raphson_7,newton_raphson_8,
    newton_raphson_9,newton_raphson_10)
  
}

probability <- function(s_target,s,volatility,delta_t_year,mean_log_variation) {
  x <- (1/(volatility*(delta_t_year)^(.5)))*(log(s_target/s)-(mean_log_variation-.5*(volatility)^(2))*delta_t_year)
  1-pnorm(x,0,1,lower.tail = TRUE)
}
#define s_target
probability_strip <- mapply(
  probability,
  350,s,volatility,delta_t_year,mean_log_variation
)

initial_volatility <- round(volatility,digits=1)

  implied_volatility_strip <- mapply(
    implied_volatility,
    initial_volatility,delta_t_year,s,strike_final,.03,tail(instrument_2[,"price"],prices_length-parameter_days_mean)
  )
  
  implied_volatility_matrix <- matrix(unlist(implied_volatility_strip),
                                      ncol=prices_length-parameter_days_mean,
                                      byrow = FALSE
  )
  
  
  k_stock <- 145.05
  s_target <- 150
  
  z_prob_itm_expiry <- (1/(volatility*(delta_t_year)^(.5)))*(log(s_target/s)-(mean_log_variation_final-.5*(volatility)^(2))*delta_t_year)
  
  prob_itm_expiry <- 1 - mapply(pnorm,z_prob_itm_expiry,0,1,lower.tail = TRUE)
  

  
  
  
  
write.xlsx(
            cbind(
              tail(instrument_2[,"effective_date"],prices_length-parameter_days_mean),
              volatility,s,mean_log_variation,delta_t_year,
              calls_prices,
              puts_prices,
              greek_delta_call,
              greek_delta_put,
              greek_gamma_strip,
              greek_theta_call,
              greek_theta_put,
              greek_vega_strip,
              greek_rho_call,
              greek_rho_put,
              probability_strip
              ),'D:/a/inversion/run_results/202001/calls_aapl_500_20200320_3.xlsx', sheetName="results",col.names=TRUE, row.names=TRUE, append=FALSE, showNA=FALSE)