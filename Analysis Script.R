#https://www.cdc.gov/nchs/nhis/2019nhis.htm
#https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/NHIS/2019/adult-codebook.pdf
#https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/NHIS/2019/adult-summary.pdf



####################################################################
####################################################################
## Outline
##
## 1) Packages and Setup
## 2) Backward Deletion
## 3) Repeat the Model Selection Protocol
## 4) Simulated Data Model Selection
## 5) Real Data Preparation
## 6) Real Data Model Selection



####################################################################
## 1) Packages and Setup
library("dplyr")
library("stringr")

"%!in%" <- function(x,y)!('%in%'(x,y))

one.mult <- function(range){
  # randomly generate a selection from a string of values
  # simulating a randomized multinomial distribution
  # where all values are treated as though equally likely
  result <- sample(range, size = 1, replace = TRUE)
  return(result)
}


modulus <- function(y, mu){
  # calculate the modulus ||a||^2
  result = (y - mu)^2 %>% sum()
  return(result)
}


add_error <- function(data, error_sd){
  # generate the y values with random error
  y <- c()
  for(i in 1:nrow(data)){
    error = rnorm(1, mean = 0, sd = error_sd)
    y[i] = data$y_star[i] + error
  }
  
  cbind(y, data)
}


tallyVar <- function(repeateResult){
  # separates the back deleting results in the repeated experiments
  # and counts the number of times a variable shows up
  repeateResult[[1]] %>% str_split(., " ~ ") %>% 
    unlist() %>% .[. != "y"]  %>% str_split(., " \\+ ") %>% 
    unlist %>% table()
}


a = test_data$y - test_data$y_star
(1/(length(a) - 1))*sum((a - mean(a))^2)

(1/(length(test_data$y) - 1))*sum((test_data$y - mean(test_data$y))^2)

####################################################################
## 2) Backward Deletion

backdelete <- function(data, me_method = c("Mallows", "LilBoot", "ConditionalBoot"), stop_threshold){
  candidate_variables = names(data)[-c(1:2)]
  candidate_remove = candidate_variables
  
  
  remove_var = NULL
  k = 1
  while(k <= length(candidate_variables)){
    # setup the variables to include in the OLS equation
    J = length(candidate_remove)
    
    
    # ME_hat for the full model
    m_eq = str_c("y", str_c(candidate_variables, collapse = " + "), sep = " ~ ")
    m_model = lm(m_eq, data)
    
    m_rss = modulus(data$y, predict(m_model)) # RSS general form
    
    
    if(k == 1){
      # S^2 estimate for sigma^2
      isolate_error = abs(data$y - predict(m_model))
      s_squared = (1/(length(isolate_error) - 1))*
                  sum((isolate_error - mean(isolate_error))^2)
      
    }
    
    
    # ME_hat before current level of subsetting
    pre_eq = str_c("y", str_c(candidate_remove, collapse = " + "), sep = " ~ ")
    pre_model = lm(pre_eq, data)
    
    pre_rss = modulus(data$y, predict(pre_model)) # RSS general form
    
    
    if(me_method == "Mallows"){
      # Mallow's Cp
      pre_candidate_me_hat = pre_rss + (2*J - nrow(data))*s_squared
      
    } else if(me_method == "LilBoot"){
      # Little Bootstrap assuming J was not data selected
      pre_candidate_me_hat = pre_rss - m_rss + ncol(data)*s_squared + 2*s_squared*(ncol(data) - J)
      
    }
    
    
    # backwards stepwise selection
    candidate_me_hat = c()
    for(i in 1:J){
      J_sub = J - 1
      
      if(me_method == "ConditionalBoot"){
        
        rss_all = c()
        me_all = c()
        for(k in 1:5){
          data$y_tild <- predict(m_model) + rnorm(nrow(data), mean = 0, sd = s_squared)
          
          eq = str_c("y_tild", str_c(candidate_remove[-i], collapse = " + "), sep = " ~ ")
          model = lm(eq, data)
          
          # ME_hat method
          me_all[k] <- modulus(predict(m_model), predict(model))
          
        }
        
      } else{
        eq = str_c("y", str_c(candidate_remove[-i], collapse = " + "), sep = " ~ ")
        model = lm(eq, data)
      
        # RSS general form
        rss = modulus(data$y, predict(model))
        
        # True ME
        true_me = modulus(predict(m_model), predict(model))
        
      }
      
      # ME_hat method
      if(me_method == "Mallows"){
        # Mallow's Cp 
        candidate_me_hat[i] = rss + (2*J_sub - nrow(data))*s_squared
        
      } else if(me_method == "LilBoot"){
        # Little Bootstrap assuming J was not data selected
        candidate_me_hat[i] = rss - m_rss + ncol(data)*s_squared + 2*s_squared*(ncol(data) - J_sub)
      
      } else if(me_method == "ConditionalBoot"){
        # Conditional Bootstrap
        me_var = var(me_all)
        candidate_me_hat[i] = mean(me_all)
        
      }
      
    }
    
    
    # check model improves
    if(me_method == "ConditionalBoot"){
      improvement = me_var < stop_threshold
      
    } else{
      improvement = abs(true_me - min(candidate_me_hat))/true_me < stop_threshold
      
    }
    
    
    if(improvement == TRUE){
      selected_var = candidate_remove
      return(selected_var)
      
      # exit function
      break
      
    }
    
    
    # minimized ME_hat to remove
    remove_var = which(candidate_me_hat == min(candidate_me_hat))[1]
    candidate_remove <- candidate_remove[-remove_var]
    
    # increase the while loop counter
    k <- k + 1
    
  }
}

# add N(0, 2) error to the y values
# Included only to test below, but excluded for the repeatBackDelete() function
# test_data = add_error(result60, 2)

# backdelete(test_data, me_method = "Mallows", 0.02)
# backdelete(test_data, me_method = "LilBoot", 0.02)
# backdelete(test_data, me_method = "ConditionalBoot", 0.75)



####################################################################
## 3) Repeat the Model Selection Protocol

repeatBackDelete <- function(dataset, method = c("Mallows", "LilBoot", "ConditionalBoot"), 
                             threshold, replicate, add_y = TRUE){
  variations = list()
  result = list()
  for(i in 1:replicate){
    
    if(add_y == TRUE){
      # add N(0, 2) error to the y values
      data = add_error(dataset, 2)
    } else if(add_y == FALSE){
      data = dataset
    }
    
    # fit with the actual model variables
    mu_star = lm(y ~ x01 + x02 + x03 + x04 + x07 + x20 + x29 + x31 + x36, data)
    star_prediction = predict(mu_star)
    
    # backward deletion to select the best model
    min_model = backdelete(data, me_method = method, threshold)
    
    eq = str_c("y", str_c(min_model, collapse = " + "), sep = " ~ ")
    model_result = lm(eq, data)
    
    variations[[i]] <- eq
    result[[i]] <- data.frame("RSS" = modulus(data$y, predict(model_result)),
                              "Actual_ME" = modulus(predict(model_result), star_prediction))
  }
  models <- do.call(rbind, variations)
  table <- do.call(rbind, result)
  
  list(models, data.frame("RSS_var" = round(apply(table[1], 2, mean), digits = 2),
             "RSS_mean" = round(apply(table[1], 2, var), digits = 2),
             "Actual_ME_var" = round(apply(table[2], 2, mean), digits = 2),
             "Actual_ME_mean" = round(apply(table[2], 2, var), digits = 2)) )
}



####################################################################
## 4) Simulated Data Model Selection
# y_star = 3*x_vec$x01 - 4.59*x_vec$x02 + 2.76*x_vec$x03 + x_vec$x04 - x_vec$x07 +
#          10*x_vec$x20 - 3.37*x_vec$x29 + 2*x_vec$x31 + x_vec$x36

fixed_x <- function(samples){
  ideal_model <- NULL
  for(i in 1:samples){
    # variables with no covariance
    x_ind <- data.frame(
      "x01" = rnorm(1, mean = 0, sd = 1),
      "x02" = rpois(1, lambda = 2),
      "x03" = one.mult(1:6),
      "x05" = rexp(1, rate = 2),
      "x06" = rexp(1, rate = 1),
      "x07" = rnorm(1, mean = 4, sd = 1.3),
      "x08" = one.mult(1:6),
      "x09" = one.mult(c(1, 2, 7, 8, 9)),
      "x10" = one.mult(1:2),
      "x12" = rnorm(1, mean = 0, sd = 2),
      "x13" = rnorm(1, mean = 0, sd = 2.5),
      "x15" = rnorm(1, mean = 2, sd = 3),
      "x16" = rnorm(1, mean = 2, sd = 1),
      "x17" = rexp(1, rate = 0.8),
      "x18" = one.mult(c(1, 5, 8, 9)),
      "x19" = rpois(1, lambda = 0.8),
      "x20" = rpois(1, lambda = 1),
      "x21" = rpois(1, lambda = 1.5),
      "x22" = rnorm(1, mean = 2.5, sd = 1),
      "x23" = rnorm(1, mean = 2.5, sd = 2),
      "x24" = rexp(1, rate = 0.6),
      "x25" = rexp(1, rate = 1.2),
      "x26" = rexp(1, rate = 1.4),
      "x27" = rexp(1, rate = 2.2),
      "x28" = rexp(1, rate = 2.4),
      "x29" = rnorm(1, mean = 4, sd = 0.5),
      "x30" = rnorm(1, mean = 5, sd = 0.75),
      "x32" = rnorm(1, mean = 5, sd = 1),
      "x33" = rnorm(1, mean = 5, sd = 1.5),
      "x34" = rnorm(1, mean = 5, sd = 2),
      "x35" = rnorm(1, mean = 5, sd = 3),
      "x36" = one.mult(c(1, 1, 1, 2, 3, 4, 4, 4, 5)),
      "x37" = rpois(1, lambda = 3),
      "x38" = rnorm(1, mean = 3, sd = 2),
      "x39" = rnorm(1, mean = 3, sd = 1),
      "x40" = rnorm(1, mean = 3, sd = 0.5)
    )
    
    # variables with covariance
    
    # generates co-linear variable to one of the multinomial variables
    if(x_ind$x08 == 1 | x_ind$x08 == 3 | x_ind$x08 == 5){
      colinear_x08 = 1
    } else if(x_ind$x08 == 2 | x_ind$x08 == 4 | x_ind$x08 == 6){
      colinear_x08 = 0
    }
    
    x_covar <- data.frame(
      "x04" = 5.34*x_ind$x05 + 2*x_ind$x01,
      "x11" = colinear_x08,
      "x14" = 5.34*x_ind$x06 + 2*x_ind$x12,
      "x31" = 5.34*x_ind$x24 + 2*x_ind$x29 + 4.56*x_ind$x20
    )
    
    # combined vector
    x_vec <- cbind(x_ind, x_covar)
    reorder = str_c("x", c("01", "02", "03", "04", "05", "06", 
                           "07", "08", "09", 10:40), sep = "")
    
    x_vec <- x_vec[order(match(names(x_vec), reorder))]
    
    # real model
    y_star = 3*x_vec$x01 - 4.59*x_vec$x02 + 2.76*x_vec$x03 + x_vec$x04 - x_vec$x07 +
             10*x_vec$x20 - 3.37*x_vec$x29 + 2*x_vec$x31 + x_vec$x36
    
    # construct the initial x-fixed data set with the real y-values
    ideal_model[[i]] <- data.frame(y_star, x_vec)
  }
  iModel = do.call(rbind, ideal_model)
  
  iModel
}


result60 = fixed_x(60)

model_Mallows60 <- repeatBackDelete(result60, method = "Mallows", 
                                    threshold = 0.02, replicate = 20)
model_LilBoot60 <- repeatBackDelete(result60, method = "LilBoot", 
                                    threshold = 0.02, replicate = 20)
model_CondBoot60 <- repeatBackDelete(result60, method = "ConditionalBoot", 
                                     threshold = 0.75, replicate = 20)

tallyVar(model_Mallows60)
tallyVar(model_LilBoot60)
tallyVar(model_CondBoot60)

table(model_Mallows60[[1]]) %>% sort(., decreasing = TRUE) %>% .[1:2]
table(model_LilBoot60[[1]]) %>% sort(., decreasing = TRUE) %>% .[1:2]
table(model_CondBoot60[[1]]) %>% sort(., decreasing = TRUE) %>% .[1:2]



result300 = fixed_x(300)

model_Mallows300 <- repeatBackDelete(result300, method = "Mallows",
                                     threshold = 0.02, replicate = 20)
model_LilBoot300 <- repeatBackDelete(result300, method = "LilBoot",
                                     threshold = 0.02, replicate = 20)
model_CondBoot300 <- repeatBackDelete(result300, method = "ConditionalBoot", 
                                      threshold = 0.75, replicate = 20)



####################################################################
## 5) Real Data Preparation
raw_data = read.csv(file.choose())

# remove the columns with 50 or more NA's in the whole array
less_50_NAs = c(1:534)[apply(raw_data, 2, function(x) sum(is.na(x))) < 50]
raw_data = raw_data[, less_50_NAs]


# select specific variables
keepVar = c("DIBEV_A", "PREDIB_A", "COPDEV_A", "ARTHEV_A", "DEMENEV_A", 
            "WEIGHTLBTC_A", "BMICAT_A", "ANXEV_A", "DEPEV_A", "HEIGHTTC_A",
            "AGEP_A", "RATCAT_A")

# randomly select remaining variables
randVar = names(raw_data)[names(raw_data) %!in% 
                          keepVar][sample(1:158, 29, replace = FALSE)]

data <- raw_data[, c(keepVar, randVar)]

dataNames <- data
names(dataNames) <- c("y", str_c("x", c("01", "02", "03", "04", "05", "06", 
                      "07", "08", "09", 10:40), sep = ""))

dataNames$y_star <- rep(1, nrow(dataNames))


reorder = c("y", "y_star", str_c("x", c("01", "02", "03", "04", "05", "06", 
                    "07", "08", "09", 10:40), sep = ""))
dataNames <- dataNames[, order(match(names(dataNames), reorder))]



####################################################################
## 6) Real Data Model Selection
# y = DIBEV_A and y_star = random, to be ignored

data60 = dataNames[sample(1:158, 60, replace = FALSE), ]

model_Mallows_real <- repeatBackDelete(data60, method = "Mallows", threshold = 0.02, 5, add_y = FALSE)
model_LilBoot_real <- repeatBackDelete(data60, method = "LilBoot", threshold = 0.3, 5, add_y = FALSE)
model_CondBoot_real <- repeatBackDelete(data60, method = "ConditionalBoot", threshold = 0.02, 5, add_y = FALSE)

tallyVar(model_Mallows_real)
tallyVar(model_LilBoot_real)
tallyVar(model_CondBoot_real)



model_Mallows_realAll <- repeatBackDelete(dataNames, method = "Mallows", threshold = 0.02, 5, add_y = FALSE)
model_LilBoot_realAll <- repeatBackDelete(dataNames, method = "LilBoot", threshold = 0.3, 5, add_y = FALSE)
model_CondBoot_realAll <- repeatBackDelete(dataNames, method = "ConditionalBoot", threshold = 0.02, 5, add_y = FALSE)

tallyVar(model_Mallows_realAll)
tallyVar(model_LilBoot_realAll)
tallyVar(model_CondBoot_realAll)







