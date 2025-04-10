library(RMySQL)
#library(xlsx)
start_time <- Sys.time()
instrument_code_init<-"'AAPL210716C00095000'"
option_database<-dbConnect(MySQL(),user='root',password='grha1j9d',dbname='options',host='localhost')

instrument<-dbGetQuery(option_database,
			paste("select * from options.extract_instrument_analysis_3 where instrument_code_split=",
				instrument_code_init,
				"order by effective_date",collapse=NULL
				)
)

future_projection_days<-dbGetQuery(option_database,
				paste("select fecha_id from options.master_operaciones where fecha_id>(select max(effective_date) from options.extract_instrument_analysis_3 where instrument_code_split=",
					instrument_code_init,
					")",
					"and fecha_id<=(select min(expiration_date) from options.extract_instrument_analysis_3 where instrument_code_split=",
					instrument_code_init,
					") and cboe=1"
					)
)

stock_price<-instrument[,"root_price"]

prices_length<-length(stock_price)

parameter_days_mean<-20

stock_prices<-lapply((parameter_days_mean+1):prices_length,
			function(i) {stock_price[(i-parameter_days_mean):i]}
)

stock_prices_2<-matrix(unlist(stock_prices),
			  ncol=prices_length-parameter_days_mean,
			  byrow=FALSE
)

stock_prices_lag<-stock_prices_2[1:(nrow(stock_prices_2)-1),,drop=FALSE]

stock_prices_2<-stock_prices_2[-1,,drop=FALSE]

log_variation<-log(stock_prices_2*(1/stock_prices_lag))

mean_log_variation_suma <- apply(log_variation,
				  2,
				  sum
)

mean_log_variation_length <- apply(log_variation,
				   2,
				   length
)

mean_log_variation <- mean_log_variation_suma/mean_log_variation_length

mean_matrix <- matrix(unlist(lapply(1:(prices_length-parameter_days_mean),
				    function(i) {rep(mean_log_variation[i],20)})
				),
				ncol=(prices_length-parameter_days_mean),
				byrow = FALSE
)

mean_substraction <- log_variation-mean_matrix

mean_sq <- mean_substraction^2

sum_sq_1 <- apply(mean_sq,
		  2,
		  sum
)

sum_sq_2 <- apply(mean_sq,
		  2,
		  length
)

sum_sq_3 <- sum_sq_1/sum_sq_2

volatility_actual <- ((sum_sq_3)^(.5))*(252^(.5))

s_actual <- tail(stock_price,
		prices_length-parameter_days_mean
)

remaining_days <- length(unlist(future_projection_days))

strike_final <- max(ifelse(
                            is.na(instrument[,"strike"]),
                    0,
                    instrument[,"strike"]
                    )
)

expiration_date <- max(ifelse(
                                is.na(instrument[,"expiration_date"]),
                       0,
                       instrument[,"expiration_date"]*1
                       )
)

delta_t_year <- unlist(lapply(1:(prices_length-parameter_days_mean+remaining_days),
				function(i) {sum(rep(1,prices_length-parameter_days_mean+remaining_days-i))}
			)
)*(1/252)

remaining_time_actual<-head(delta_t_year,prices_length-parameter_days_mean)

remaining_time_to_expiration<-tail(delta_t_year,remaining_days)

black_scholes <- function(S,sigma,delta_t,K,r) {
  d1 <- (1/(sigma*(delta_t)^(.5)))*(log(S/K)+(r+.5*(sigma^2))*delta_t)
  d2 <- d1 - sigma*(delta_t^(.5))
    black_scholes_call<-S*pnorm(d1,0,1,lower.tail = TRUE)-K*pnorm(d2,0,1,lower.tail = TRUE)*exp(-r*delta_t)
    black_scholes_put<-K*pnorm(-d2,0,1,lower.tail = TRUE)*exp(-r*delta_t)-S*pnorm(-d1,0,1,lower.tail = TRUE)
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

black_scholes_strip <- mapply(black_scholes,
                              s_actual,
			      volatility_actual,
			      remaining_time_actual,
			      strike_final,
			      .03
)

black_scholes_prices_matrix <- matrix(unlist(black_scholes_strip),
                                      ncol=(prices_length-parameter_days_mean)
)

calls_prices <- black_scholes_prices_matrix[1,]

puts_prices <- black_scholes_prices_matrix[2,]


greek_delta_strip <- mapply(greek_delta,
                            s_actual,
			    volatility_actual,
			    remaining_time_actual,
			    strike_final,
			    .03
)

end_time <- Sys.time()
total_execution_time <- end_time - start_time
# 20230121 llegue hasta aca testeandolo en R
