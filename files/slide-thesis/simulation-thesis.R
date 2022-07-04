# Data
pacman::p_load(tidyverse)
pacman::p_load(plotly)

# models

# --------

fitLocalModels <- function(kmeans_BD,
                           train_response,
                           model = "lm",
                           formula = NULL){
  start_time <- Sys.time()
  X_train <- kmeans_BD$train_data$X_train
  y_train <- train_response[!(kmeans_BD$train_data$id_remain)]
  X_remain <- kmeans_BD$train_data$X_remain
  nr_remain <- nrow(X_remain)
  y_remain <- NULL
  if(nr_remain != 0){
    y_remain <- train_response[kmeans_BD$train_data$id_remain]
  }
  pacman::p_load(tree)
  pacman::p_load(randomForest)
  model_ <- ifelse(model == "tree", tree::tree, model)
  K <- nrow(kmeans_BD$centroids)
  if (is.null(formula)){
    form <- formula(target ~ .)
  }
  else{
    form <- update(formula, target ~ .)
  }
  data_ <- bind_cols(X_train, "target":= y_train)
  fit_lookup <- list(lm = "fitted.values",
                     rf = "predicted")
  if(is.character(model_)){
    model_lookup <- list(lm = lm,
                         rf = randomForest::randomForest)
    mod <- map(.x = 1:K, 
               .f = ~ model_lookup[[model_]](formula = form, 
                                             data = data_[kmeans_BD$clusters == .x, ]))
  } else{
    mod <- map(.x = 1:K, 
               .f = ~ model_(formula = form, 
                             data = data_[kmeans_BD$clusters == .x,]))
  }
  pred0 <- NULL
  if(nr_remain != 0){
    pred0 <- vector(mode = "numeric", 
                    length = length(y_remain))
    clus <- findClosestCentroid(x. = X_remain,
                                centroids. = kmeans_BD$centroids,
                                div = kmeans_BD$parameters$div,
                                deg = kmeans_BD$parameters$deg)
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(mod[[i_]],
                                   as.data.frame(X_remain[clus == i_,]))
    }
  }
  time_taken <- Sys.time() - start_time
  return(list(
    local_models = mod,
    kmeans_BD = kmeans_BD,
    data_remain = list(fit = pred0,
                       response = y_remain),
    running_time = time_taken
  ))
}


localPredict <- function(localModels,
                         newData){
  kmean_BD <- localModels$kmeans_BD
  K <- nrow(kmean_BD$centroids)
  newData_ <- newData
  if(!(is.null(kmean_BD$parameters$center_))){
    newData_ <- scale(newData,
                      center = kmean_BD$parameters$center_,
                      scale = kmean_BD$parameters$scale_)
    id0 <- (newData_ <= 0)
    if(sum(id0) > 0){
      min_ <- min(newData_[id0])
      newData_[id0] <- runif(sum(id0), min(1e-3, min_/10), min_)
    }
  }
  clus <- findClosestCentroid(x. = newData_,
                              centroids. = kmean_BD$centroids,
                              div = kmean_BD$parameters$div,
                              deg = kmean_BD$parameters$deg)
  pred0 <- vector(mode = "numeric", length = nrow(newData_))
  for(i_ in 1:K){
    pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                 as.data.frame(newData_[clus == i_,]))
  }
  pred0 <- as_tibble(pred0)
  names(pred0) <- ifelse(kmean_BD$parameters$div == "polynomial",
                         paste0("polynomial", kmean_BD$parameters$deg),
                         kmean_BD$parameters$div)
  return(pred0)
}

dataSimulation = function(m, sigma, prop.test = 0.2, distrib, model, figure = FALSE){
  start.time <- Sys.time()
  DIM3 = F
  LOG = F
  M = round((1 + prop.test) * m)
  sigmoid = function(x){
    return(1/(1+exp(-x)))
  }
  if((distrib == "Exp") || (distrib == "exp") ||(distrib == "exponential") || 
     (distrib == "Exponential")){
    x1 = c(rexp(M, 0.05) + 9, rexp(M, 0.5), rexp(M,0.1) + 8)
    x2 = c(rexp(M, 0.5), rexp(M, 0.05), rexp(M,0.1) + 8)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "Pois") || (distrib == "pois") ||(distrib == "poisson") || 
     (distrib == "poisson")){
    x1 = c(rpois(M, 3), rpois(M,10), rpois(M,13) + 2)
    x2 = c(rpois(M, 11), rpois(M,2), rpois(M,12) + 2)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "geom") || (distrib == "Geom") ||(distrib == "Geometric") || 
     (distrib == "geometric")){
    x1 = c(rgeom(M, 0.07)+9, rgeom(M, 0.55),rgeom(M, 0.15)+7)
    x2 = c(rgeom(M, 0.35), rgeom(M, 0.07)+17,rgeom(M, 0.15)+12)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "gaussian2") || (distrib == "Gaussian2") ||(distrib == "Normal2") || 
     (distrib == "normal2")){
    x1 = c(rnorm(M, 4, 1), rnorm(M, 22, 2),rnorm(M, 10, 2))
    x2 = c(rnorm(M, 12, 1), rnorm(M, 9, 1),rnorm(M, 5, 2))
    x = cbind(x1, x2)
    x = x + abs(min(x)) + 1
  }
  if((distrib == "gaussian3") || (distrib == "Gaussian3") ||(distrib == "Normal3") || 
     (distrib == "normal3")){
    x1 = c(rnorm(M, 6, 1), rnorm(M, 5, 2),rnorm(M, 8, 1))
    x2 = c(rnorm(M, 14, 2), rnorm(M, 10, 1),rnorm(M, 6, 1))
    x3 = c(rnorm(M, 6, 1), rnorm(M, 15, 2),rnorm(M, 14, 2))
    x = cbind(x1,x2,x3) 
    x = x + abs(min(x)) + 5
    DIM3 = T
  }
  if((model == "logit") || (model == "logistic") || (model == "Logit") || 
     (model == "Logistic")){
    LOG = T
  }
  id.train = c(1:m,(M+1):(M+m),(2*M+1):(2*M+m))
  x.train = x[id.train,]
  x.test = x[-id.train,]
  df = as.matrix(cbind(x0 = 1, x))
  p = ncol(x)
  if(DIM3){
    b1 = c(-10,3,7)
    b2 = c(7,5,-12)              # We fixed the coeffiencets of the 3D case
    b3 = c(6,-11,10)
    beta = c(b1,b2,b3)
  }
  else{
    b1 = c(-8,3)
    b2 = c(-6,-5)                 # We fixed the coeffiencets of the 2D case
    b3 = c(5,-7)
    beta = c(b1,b2,b3)
  }
  
  if(LOG){
    beta1 = c(-colMeans(x.train[1:m,]) %*% beta[1:p], beta[1:p])
    beta2 = c(-colMeans(x.train[(m+1):(2*m),]) %*% beta[(p+1):(2*p)], 
              beta[(p+1):(2*p)])
    beta3 = c(-colMeans(x.train[(2*m+1):(3*m),]) %*% beta[(2*p+1):(3*p)], 
              beta[(2*p+1):(3*p)])
    y.1 = df[1:M,] %*% beta1
    y.2 = df[(M+1):(2*M),] %*% beta2 
    y.3 = df[(2*M+1):(3*M),] %*% beta3
    z = sigmoid(c(y.1,y.2,y.3)+ rnorm(3*M,0,sigma))
    target = factor(as.numeric(z > 0.5), levels = c(1,0))
  }
  else{
    beta1 = c(-15,beta[1:p])
    beta2 = c(25,beta[(p+1):(2*p)])
    beta3 = c(-10,beta[(2*p+1):(3*p)])
    y.1 = df[1:M,] %*% beta1
    y.2 = df[(M+1):(2*M),] %*% beta2 
    y.3 = df[(2*M+1):(3*M),] %*% beta3
    target = c(y.1, y.2,y.3) + rnorm(3*M,0,sigma)
  }
  cluster = c(rep(1,M), rep(2,M), rep(3,M))
  if(LOG){
    y.train = target[id.train]
    y.test = target[-id.train]
    col.train = as.numeric(y.train) + 2
  }
  else{
    y.train = target[id.train]
    y.test = target[-id.train]
    col.train = c(rep(2,m), rep(3,m), rep(4,m))
  }
  if(figure){
    if(DIM3){
      par(mfrow=c(1,1))
      scatterplot3d::scatterplot3d(x.train, color = cluster[id.train] + 1, 
                                   main = "Simulated Predictor", xlab = "x1", 
                                   ylab = "x2", zlab = "x3")
    }
    else{
      par(mfrow=c(1,2))
      if(LOG){
        plot(x.train, col=cluster[id.train] + 1, main = "Predictors", xlab = "x1", ylab = "x2")
        scatterplot3d::scatterplot3d(x=x.train[,1], y=x.train[,2],z=y.train, color = as.numeric(y.train) + 5, main = "Completed Data", xlab = "x1", ylab = "x2", zlab = "y")
      }
      else{
        plot(x.train[1:m,], col=2, main = "Predictors", xlim = c(min(x.train[,1])-5, max(x.train[,1])+5), ylim = c(min(x.train[,2])-5, max(x.train[,2])+5), pch = 1, cex.main = 2)
        points(x.train[(m+1):(2*m),], col=3, pch=17)
        points(x.train[(2*m+1):(3*m),], col=4, pch=8)
        scatterplot3d::scatterplot3d(x=c(x.train[,1],x.train[,1]), y=c(x.train[,2],x.train[,2]),z=c(y.train + max(y.train), rep(min(y.train),length(y.train))), color = rep(cluster[id.train] + 1,2), main = "Complete Data", xlab = "x1", ylab = "x2", zlab = "y",  cex.main = 2, angle = 20)
      }
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  return(list(Training.Data = as.data.frame(cbind(x.train,y.train)), 
              Testing.Data = as.data.frame(cbind(x.test,y.test)), 
              BETA = cbind(beta1, beta2, beta3), 
              CLUSTERS = cluster, Running.Time = time.taken))
}

# =============


x1 <- c(runif(100,-1,0), rexp(100, 2), rnorm(100,-0.65,0.3))
x2 <- c(runif(100,-2,-0.1), rexp(100, 1), rnorm(100,1,0.4))
y <- rep(c(-0.1,-0.2, 0.3), each = 100) * x1 + rep(c(-0.2,0.1,-0.1), each = 100) * x2 + rnorm(300, 0, .05) - rep(c(0.7,.2,-.2), each = 100) + 1
col1 <- rep("group3", 300)
col1[3*x2 < x1] <- "group1"
col1[(3*x2 >= x1) & ( 2*x1 < -x2)] <- "group2"
col2 <- rep("group3", 300)
col2[1:100] <- "group1"
col2[x1^2 + x2^2 < 0.5] <- "group1"
col2[( x1 > 0.4) & (col2 != "group1")] <- "group2"
df0 <- tibble(x1 = rep(x1, 2), 
              x2 = rep(x2, 2),
              y = c(y, rep(0, 300)),
             color0 = rep(rep(c("group1", "group2", "group3"), each = 100), 2))
fig0 <- plotly::plot_ly(df0,
                        x = ~x1,
                        y = ~x2,
                        z = ~y,
                        size = .5) %>%
  plotly::add_markers(color = ~color0)
fig0

# clustering

df <- tibble(x1 = x1,
             x2 = x2)
devtools::source_url("https://raw.githubusercontent.com/hassothea/KFC-Procedure/master/kmeanBD.R")

div_list <- c("euclidean",
              "gkl",
              "logistic",
              "itakura",
              "exponential",
              "polynomial")
min_ <-  apply(df, 2, min) - 0.05
scale_ <-  apply(df, 2, max) - min_ + 0.05
df_scale <- scale(df,
                  center = min_,
                  scale = scale_)
clusters <- div_list %>% map(.f = ~ kmeansBD(train_input = df_scale,
                                 K = 3,
                                 deg = 4,
                                 div = .))
names(clusters) <- div_list

df0 %>%
  mutate(col1 = rep(col1, 2)) %>%
  mutate(col2 = rep(col2, 2)) %>%
  mutate(col_eucl = rep(as_factor(clusters$euclidean$clusters), 2)) %>%
  mutate(col_gkl = rep(as_factor(clusters$gkl$clusters), 2)) %>%
  mutate(col_log = rep(as_factor(clusters$logistic$clusters), 2)) %>%
  mutate(col_ita = rep(as_factor(clusters$itakura$clusters), 2)) %>%
  mutate(col_poly4 = rep(as_factor(clusters$polynomial$clusters), 2)) %>%
  mutate(col_exp = rep(as_factor(clusters$exponential$clusters), 2)) -> df1

fig <- plotly::plot_ly(df1,
                       x = ~x1,
                       y = ~x2,
                       z = ~y,
                       size = .5)
fig %>%
  add_markers(color = ~col_eucl)
fig %>%
  add_markers(color = ~col_gkl)
fig %>%
  add_markers(color = ~col_log)
fig %>%
  add_markers(color = ~col_ita)
fig %>%
  add_markers(color = ~col_exp)


df2d <- tibble("x1":= x1, 
               "x2":= x2,
               "y":=y,
               "col0":= rep(paste0("group", 1:3), each = 100),
               "col1":= col1,
               "col2":=col2,
               "euc":= factor(clusters$euclidean$clusters, labels = paste0("group", 1:3)),
               "gkl":= factor(clusters$gkl$clusters, labels = paste0("group", 1:3)),
               "ita":= factor(clusters$itakura$clusters, labels =  paste0("group", 1:3)),
               "log":= factor(clusters$logistic$clusters, labels = paste0("group", 1:3)),
               "poly4":= factor(clusters$polynomial$clusters, labels = paste0("group", 1:3)),
               "exp":= factor(clusters$exponential$clusters, labels =  paste0("group", 1:3)))

fig2d <- plot_ly(data = df2d,
                 x = ~x1,
                 y = ~x2,
                 size = .5,
                 showlegend = FALSE)
fig2d %>%
  add_markers(color = ~ ita) -> fig2d1
fig2d %>%
  add_markers(color = ~ col2) -> fig2d2
fig2d %>%
  add_markers(color = ~ col0) -> fig2d3
subplot(fig2d1, fig2d2, fig2d3, nrows = 1)

# fitting models

models <- clusters %>%
  map(.f = ~ fitLocalModels(.x, y))

preds <- models %>%
  map(.f = ~ localPredict(.x, newData = df_scale))

df2d_ <- df2d %>%
  mutate(pred_euc = preds$euclidean$euclidean) %>%
  mutate(pred_gkl = preds$gkl$gkl) %>%
  mutate(pred_log = preds$logistic$logistic) %>%
  mutate(pred_ita = preds$itakura$itakura) %>%
  mutate(pred_exp = preds$exponential$exponential) %>%
  mutate(pred_poly4 = preds$polynomial$polynomial4)


# with predictions
df3d <- df1 %>%
  mutate(pred_euc = c(preds$euclidean$euclidean, rep(0, 300))) %>%
  mutate(pred_gkl = c(preds$gkl$gkl, rep(0, 300))) %>%
  mutate(pred_log = c(preds$logistic$logistic, rep(0, 300))) %>%
  mutate(pred_ita = c(preds$itakura$itakura, rep(0, 300))) %>%
  mutate(pred_exp = c(preds$exponential$exponential, rep(0, 300))) %>%
  mutate(pred_poly4 = c(preds$polynomial$polynomial4, rep(0, 300)))



# -------------------------------- clustering -----------------------------


df_norm <- dataSimulation(200, 5, 0.3, distrib = "normal2", model = "reg" ,figure = F)
df_pois <- dataSimulation(200, 5, 0.3, distrib = "pois", model = "reg" ,figure = F)
df_exp <- dataSimulation(200, 5, 0.3, distrib = "exp", model = "reg" ,figure = F)

list_X <- list(norm = df_norm$Training.Data[,1:2],
               pois = df_pois$Training.Data[,1:2],
               exp = df_exp$Training.Data[,1:2])

div_list <- c("euclidean",
              "gkl",
              "polynomial",
              "itakura")

clust <- list()

for(x_ in c("norm", "pois", "exp")){
  #min_ <-  apply(list_X[[x_]], 2, min) - 0.05
  #scale_ <-  apply(list_X[[x_]], 2, max) - min_ + 0.05
  # df_scale <- scale(list_X[[x_]],
  #                   center = min_,
  #                   scale = scale_)
  clust[[x_]] <- div_list %>%
    map(.f = ~kmeansBD(train_input = list_X[[x_]], 
                       K = 3,
                       deg = 5,
                       div = .x)$clusters)
  names(clust[[x_]]) <- div_list
}

color0 <- as_factor(rep(1:3, each = 200))
df0 <- tibble(x1 = list_X[["exp"]]$x1, 
              x2 = list_X[["exp"]]$x2,
              color0 = color0,
              color_euc = as_factor(clust$exp$euclidean),
              color_gkl = as_factor(clust$exp$gkl),
              color_poly5 = as_factor(clust$exp$polynomial),
              color_ita = as_factor(clust$exp$itakura))
fig0 <- plot_ly(df0,
                x = ~x1,
                y = ~x2,
                showlegend = FALSE)

fig1 <- fig0 %>%
  add_markers(color = ~color_poly5,
              size = I(150)) %>%
  layout(title = list(text = "Euclidean divergence",
                      y = 0.925))
fig1


x <- seq(0.1, 5, length.out = 200)
y <- -exp(-x)*log(x^20) + 4
x_m <- 1.763
y_m <- -exp(-x_m)*log(x_m^20) + 4
library(magrittr)

df_ <- tibble::tibble(error = y,
                       parameter = x)
plotly::plot_ly(df_,
                x = ~parameter,
                y = ~error,
                showlegend = FALSE) %>%
  plotly::add_lines() %>%
  plotly::add_lines(x = x_m,
                    colors = "red") %>%
  plotly::add_markers(x = x_m,
                      y = y_m) %>%
  plotly::layout(scene = list(aspectration=list(x=1,y=1)))
  

# -------------------------------------------------------------------------

# saving data
write.table(as.data.frame(df1), 
      "D:/PhD-Works/Conferences/SFdS2022/slide-sfds22/slides/data/df_full1.txt",
      col.names = TRUE,
      row.names = FALSE)

write.table(as.data.frame(df2d), 
            "D:/PhD-Works/Conferences/SFdS2022/slide-sfds22/slides/data/df2d.txt",
            col.names = TRUE,
            row.names = FALSE)

write.table(as.data.frame(df2d_),
            "D:/PhD-Works/Conferences/SFdS2022/slide-sfds22/slides/data/df2d_predict.txt",
            col.names = TRUE,
            row.names = FALSE)

write.table(as.data.frame(df3d), 
            "D:/PhD-Works/Conferences/SFdS2022/slide-sfds22/slides/data/df3d_predict.txt",
            col.names = TRUE,
            row.names = FALSE)

write.table(as.data.frame(df0), 
            "D:/PhD-Works/Conferences/SFdS2022/slide-sfds22/slides/data/df_exp.txt",
            col.names = TRUE,
            row.names = FALSE)
