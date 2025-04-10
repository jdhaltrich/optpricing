library(RMySQL)

# inputs

start_time <- Sys.time()

instrument_code_init <- "'AAPL210917C00062500'"
instrument_code_export <- "AAPL210917C00062500" 
entry_price_cap <- .3
min_required_performance_pc <- 4

# calculation

st_time_db <- Sys.time()

option_database <- dbConnect(MySQL(),user='root',password='grha1j9d',dbname='options',host='127.0.0.1')

st_time_db_extract <- Sys.time()

instrument <- dbGetQuery(option_database,
			paste("select * from options.IA_BASE where instrument_code_split=",
				instrument_code_init,
				collapse=NULL
				)
)

end_time_db_extract <- Sys.time()
exec_time_db_extract <- end_time_db_extract - st_time_db_extract

st_time_db_proj <- Sys.time()

future_projection_days <- dbGetQuery(option_database,
				paste("select fecha_id from options.master_operaciones where fecha_id>(select max(effective_date) from options.IA_BASE where instrument_code_split=",
					instrument_code_init,
					")",
					"and fecha_id<=(select min(expiration_date) from options.IA_BASE where instrument_code_split=",
					instrument_code_init,
					") and cboe=1"
					)
)

end_time_db_proj <- Sys.time()
exec_time_db_proj <- end_time_db_proj - st_time_db_proj

end_time_db <- Sys.time()
exec_time_db <- end_time_db - st_time_db

stock_price <- instrument[,"root_price"]

prices_length <- length(stock_price)

parameter_days_mean <- 20

length_actual_val_days <- prices_length - parameter_days_mean

stock_prices <- lapply((parameter_days_mean+1):prices_length,
			function(i) {stock_price[(i-parameter_days_mean):i]}
)

stock_prices_2 <- matrix(unlist(stock_prices),
			  ncol=length_actual_val_days,
			  byrow=FALSE
)

stock_prices_lag <- stock_prices_2[1:(nrow(stock_prices_2)-1),,drop=FALSE]

stock_prices_2 <- stock_prices_2[-1,,drop=FALSE]

log_variation <- log(stock_prices_2*(1/stock_prices_lag))

mean_log_variation_suma <- apply(log_variation,
				  2,
				  sum
)

mean_log_variation_length <- apply(log_variation,
				   2,
				   length
)

mean_log_variation <- mean_log_variation_suma/mean_log_variation_length

mean_matrix <- matrix(unlist(lapply(1:length_actual_val_days,
				    function(i) {rep(mean_log_variation[i],20)})
				),
				ncol=length_actual_val_days,
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
		length_actual_val_days
)

length_projected_val_days <- length(unlist(future_projection_days))

total_proj_period <- length_actual_val_days + length_projected_val_days

strike_final <- tail(
			ifelse(
				is.na(instrument[,"strike"]),
				0,
				instrument[,"strike"]
			),
			length_actual_val_days
)

expiration_date <- max(ifelse(
                                is.na(instrument[,"expiration_date"]),
                       0,
                       instrument[,"expiration_date"]*1
                       )
)

delta_t_year <- unlist(lapply(1:total_proj_period,
				function(i) {sum(rep(1,total_proj_period-i))}
			)
)*(1/252)

remaining_time_actual <- head(delta_t_year,length_actual_val_days)

remaining_time_to_expiration <- tail(delta_t_year,length_projected_val_days)

proj_type <- rep("Actual_Valuation",length(s_actual))

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

black_scholes_strip <- mapply(black_scholes,
                              s_actual,
			      volatility_actual,
			      remaining_time_actual,
			      strike_final,
			      .03
)

black_scholes_prices_matrix <- matrix(unlist(black_scholes_strip),
                                      ncol=length_actual_val_days
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

greek_delta_matrix <- matrix(unlist(greek_delta_strip),
                             ncol=length_actual_val_days
)

greek_delta_call <- greek_delta_matrix[1,]

greek_delta_put <- greek_delta_matrix[2,]

greek_gamma_strip <- mapply(greek_gamma,
                            s_actual,
			    volatility_actual,
			    remaining_time_actual,
			    strike_final,
			    .03
)

greek_gamma_matrix <- matrix(unlist(greek_gamma_strip),
                             ncol=length_actual_val_days
)

greek_gamma <- greek_gamma_matrix[1,]


greek_theta_strip <- mapply(greek_theta,
                            s_actual,
			    volatility_actual,
			    remaining_time_actual,
			    strike_final,
			    .03
)

greek_theta_matrix <- matrix(unlist(greek_theta_strip),
                             ncol=length_actual_val_days
)

greek_theta_call <- greek_theta_matrix[1,]

greek_theta_put <- greek_theta_matrix[2,]


greek_vega_strip <- mapply(greek_vega,
                           s_actual,
			   volatility_actual,
			   remaining_time_actual,
			   strike_final,
			   .03
)

greek_vega_matrix <- matrix(unlist(greek_vega_strip),
                             ncol=length_actual_val_days
)

greek_vega <- greek_vega_matrix[1,]


greek_rho_strip <- mapply(greek_rho,
                          s_actual,
			  volatility_actual,
			  remaining_time_actual,
			  strike_final,
			  .03
)

greek_rho_matrix <- matrix(unlist(greek_rho_strip),
                             ncol=length_actual_val_days
)

greek_rho_call <- greek_rho_matrix[1,]

greek_rho_put <- greek_rho_matrix[2,]

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

initial_volatility <- round(volatility_actual,digits=1)
							
							
implied_volatility_strip <- mapply(implied_volatility,
                                   initial_volatility,
				   remaining_time_actual,
				   s_actual,
				   strike_final,
				   .03,
				   tail(instrument[,"price"],length_actual_val_days)
)

implied_volatility_matrix <- matrix(unlist(implied_volatility_strip),
                                    ncol=length_actual_val_days
)

implied_volatility_calculation <- implied_volatility_matrix[10,]
									

probability <- function(S_target,S,volatility,delta_t,mean_log) {
  x <- (1/(volatility*(delta_t)^(.5)))*(log(S_target/S)-(mean_log-.5*(volatility)^(2))*delta_t)
  1-pnorm(x,0,1,lower.tail = TRUE)
}


probability_strip <- mapply(probability,
                            405,
			    s_actual,
			    volatility_actual,
			    remaining_time_actual,
			    mean_log_variation
)


instrument_actuals <- cbind(
			proj_type,
			tail(instrument,length_actual_val_days)
)

rm(instrument)

actuals_calc <- cbind(
		instrument_actuals,
		volatility_actual,
		s_actual,
		mean_log_variation,
		remaining_time_actual,
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
		probability_strip,
#		prob_itm_expiry,
		implied_volatility_calculation
)

fields_actuals_calc <- colnames(actuals_calc)

actuals_calc_last <- actuals_calc[length_actual_val_days,]

fields_actuals_calc_last <- colnames(actuals_calc_last)

# s_target <- 405

# z_prob_itm_expiry <- (1/(volatility_actual*(remaining_time_actual)^(.5)))*(log(s_target/s_actual)-(mean_log_variation-.5*(volatility_actual)^(2))*remaining_time_actual)

#prob_itm_expiry <- 1 - mapply(pnorm,
#                              z_prob_itm_expiry,
#			      0,
#			      1,
#			      lower.tail = TRUE
#)            # laynestaley88: previous test to calculate GBM of future projected Underlying price

opt_target_price <- entry_price_cap*(1+ min_required_performance_pc)

# Cox, Ross and Rubinstein binomial tree pricing of American Options implementation

crr_binomial_tree <- function(S,k,u,d,p,n) {
		
}

crr_binomial_tree_american <- function(S,K,u,d,r,steps) {
	p <- (r - d)/(u - d)
	C <- c()
	for(i in 0:steps) {
		C[i] <- (factorial(steps)/((factorial(i))*(factorial(steps-i))))*(p^i)*((1 - p)^(steps - i))*pmax((u^i)*(d^(steps-i))*S-K,0)
	}
	sum(C)/(r^steps)
}


###  Start of Projected Valuation

proj_type_projected <- rep("Projected_Valuation",length_projected_val_days)

projected_calc <- cbind(proj_type_projected,rep(1,length_projected_val_days)) 

#mc_s0 <- actuals_calc_last[,fields_actuals_calc_last[29]]
#mc_sigma <- actuals_calc_last[,fields_actuals_calc_last[32]]
#mc_delta_t <- actuals_calc_last[,fields_actuals_calc_last[35]]
#mc_K <- strike_final
#mc_r <- .03


# option valuation by Montecarlo simulation

# previous implementation

runs <- 100000
generate_path <- function(underlying,media_variacion,volatility_simulation,delta_t_simulation,days) {
  change <- rnorm(days, mean = 0, sd = 1)
  sample_path <- cumprod(c(underlying,
  			   exp((media_variacion-.5*(volatility_simulation^2))*delta_t_simulation+change*volatility_simulation*(delta_t_simulation^(.5)))
			 )
		)
  list(change,sample_path,sample_path[days+1])
}

# new implementation 20230301

mc_options_eu <- function(mc_sim,mc_s0,mc_sigma,mc_delta_t,mc_K,mc_r) {
	init_time <- Sys.time()
	z <- rnorm(mc_sim,mean = 0,sd = 1)
	wt <- ((mc_delta_t)^(.5))*z
	
	mc_st <- mc_s0*exp((mc_r - .5*(mc_sigma^2))*(mc_delta_t) + mc_sigma*wt)
	mc_call_payoff <- exp(-mc_r*(mc_delta_t))*pmax(mc_st-mc_K,0)
	mc_put_payoff <- exp(-mc_r*(mc_delta_t))*pmax(mc_K-mc_st,0)
	mc_call_price <- sum(mc_call_payoff)/mc_sim
	mc_put_price <- sum(mc_put_payoff)/mc_sim
	
	mc_call_sd <- (sum(((mc_call_payoff - mc_call_price)/mc_sim)^2))^(.5)
	mc_put_sd <- (sum(((mc_put_payoff - mc_put_price)/mc_sim)^2))^(.5)
	mc_call_cv <- mc_call_sd/mc_call_price
	mc_put_cv <- mc_put_sd/mc_put_price
	end_time <- Sys.time()
	exec_time <- end_time - init_time
	#list(mc_st,mc_call_payoff,mc_put_payoff,mc_call_price,mc_put_price,init_time,end_time,exec_time)
	list(mc_call_price,mc_call_sd,mc_call_cv,mc_put_price,mc_put_sd,mc_put_cv,init_time,end_time,exec_time)
}

#mc_st_p <- function(mc_sim,mc_s0,mc_sigma,mc_delta_t,mc_K,mc_r) {
#	init_time <- Sys.time()
#	z <- rnorm(mc_sim,mean = 0,sd = 1)
#	wt <- ((mc_delta_t)^(.5))*z
	
#	mc_st_int <- mc_s0*exp((mc_r - .5*(mc_sigma^2))*(mc_delta_t) + mc_sigma*wt)
#	mc_st_final <- sum(mc_st_int)/mc_sim
#	mc_call_payoff <- exp(-mc_r*(mc_delta_t))*pmax(mc_st-mc_K,0)

#	end_time <- Sys.time()
#	exec_time <- end_time - init_time
	#list(mc_call_price,mc_call_sd,mc_call_cv,mc_put_price,mc_put_sd,mc_put_cv,init_time,end_time,exec_time)
#	list(mc_st_final)
#}

#z_test <- rnorm(10000000,mean = 0,sd = 1)
#wt_test <- ((1/252)^(.5))*z_test
	
#mc_st_int_test <- mc_s0*exp((mc_r - .5*(mc_sigma^2))*(1/252) + mc_sigma*wt_test)
#mc_st_final_test <- sum(mc_st_int_test)/10000000
#mc_call_payoff_test_p <- exp(-mc_r*(1/252))*pmax(mc_st_int_test-mc_K,0)

#mc_options_eu_proj <- for(i in 1:length_projected_val_days) {}

# Merge of Actual and Porjected Valuations

complete_proj_type <- c(
			proj_type,
			proj_type_projected
)


#0.336325799326088
#143.160004
#simulacion <- replicate(100000,generate_path(88,0.5,0.84,1/252,40))

#simulacion_precios <- simulacion[2,]

#simulacion_precios_2 <- matrix(unlist(simulacion_precios),ncol=100000,byrow = FALSE)

#simulacion_precios_3 <- simulacion_precios_2[,1:1000]

#sum(simulacion_precios_2[5,])/length(simulacion_precios_2[5,])

end_time <- Sys.time()
exec_time_model <- end_time - start_time


write.csv(actuals_calc,
		paste0("/laynestaley88/inversion/modelo/R_workspace/dev/202303/output_R_model_implementation",
			"/",
			instrument_code_export,
			"_20230314.csv",
			collapse=NULL
			)
)
