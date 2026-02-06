环境与配置
rm(list = ls())
if(!require(devEMF)) install.packages("devEMF")
library(devEMF)

total_start <- proc.time() 

#核心参数
K <- 3        
p <- 10000    
c_ratio <- 20
n <- round(p / c_ratio)  
im_part <- 0.05 


plot_xlim <- c(-3, 3)
plot_ylim <- c(0, 0.7)


data_filename <- paste0("sim_data_K", K, "_p", p, "_n", n, ".RData")
plot_file_emf <- paste0("plot_K", K, "_p", p, "_n", n, ".emf")
plot_file_png <- paste0("plot_K", K, "_p", p, "_n", n, ".png")

cat("参数设定：K=", K, ", p=", p, ", n=", n, ", p/n=", c_ratio, "\n")
cat("数据文件：", data_filename, "\n")
cat("输出图片：", plot_file_emf, " & ", plot_file_png, "\n\n")

# 模拟与数据加载模块

# 检查是否存在数据文件
if (file.exists(data_filename)) {
  cat(">>> [1/4] 检测到现有数据文件，正在加载...\n")
  load(data_filename)
  cat(">>> 数据加载完成！跳过模拟步骤。\n\n")
} else {
  cat(">>> [1/4] 未检测到数据文件，开始执行模拟...\n")
  
  #构造A矩阵
  cat("   -> 构造A矩阵...") 
  A <- vector(mode = "list", length = K)
  a2 <- numeric(K)
  for (r in 1:K) {
    lambda_A <- rexp(p, rate = 1)
    u <- runif(p); u <- u / sqrt(sum(u^2))
    P <- t(apply(matrix(1:p, ncol=1), 1, function(i, u) {
      e_i <- numeric(p); e_i[i] <- 1; e_i - 2 * u[i] * u
    }, u = u))
    # 优化构造速度
    diag_p <- P * lambda_A
    A_r <- diag_p %*% t(P) 
    A[[r]] <- A_r
    a2[r] <- mean(lambda_A^2)
  }
  cat("完成\n")
  
  #构造B矩阵
  cat("   -> 构造B矩阵...")
  B <- vector(mode = "list", length = K)
  H <- matrix(0, n, K)
  for (r in 1:K) {
    lambda_B <- rexp(n, rate = 2*r)
    H[, r] <- lambda_B
    v <- runif(n); v <- v / sqrt(sum(v^2))
    Q <- t(apply(matrix(1:n, ncol=1), 1, function(i, v) {
      e_i <- numeric(n); e_i[i] <- 1; e_i - 2 * v[i] * v
    }, v = v))
    diag_q <- Q * lambda_B
    B_r <- diag_q %*% t(Q)
    B[[r]] <- B_r
  }
  cat("完成\n")
  
  #生成Z和X矩阵
  cat("   -> 生成Z和X矩阵...")
  Z <- vector(mode = "list", length = K)
  X <- matrix(0 + 0i, p, n)
  for (r in 1:K) {
    Z[[r]] <- matrix(complex(real = rnorm(p*n)/sqrt(2), imaginary = rnorm(p*n)/sqrt(2)), nrow = p, ncol = n)
    eig_A <- eigen(A[[r]], symmetric = TRUE)
    A_sqrt <- eig_A$vectors %*% (sqrt(eig_A$values) * t(eig_A$vectors))
    eig_B <- eigen(B[[r]], symmetric = TRUE)
    B_sqrt <- eig_B$vectors %*% (sqrt(eig_B$values) * t(eig_B$vectors))
    
    A_Z <- A_sqrt %*% Z[[r]]
    X <- X + A_Z %*% B_sqrt
  }
  cat("完成\n")
  
  #计算C_n和特征值
  cat("   -> 计算特征值...")
  S_p <- (t(Conj(X)) %*% X) / p
  E_Sp <- matrix(0 + 0i, n, n)
  for (r in 1:K) {
    tr_A <- sum(diag(A[[r]]))
    E_Sp <- E_Sp + B[[r]] * (tr_A / p)
  }
  C_n <- sqrt(p / n) * (S_p - E_Sp)
  eig_Cn <- eigen(Re(C_n), only.values = TRUE)$values
  eig_Cn <- eig_Cn[abs(eig_Cn) > 1e-8] 
  
  #保存数据以便下次使用
  save(list = c("eig_Cn", "a2", "H", "K", "p", "n", "c_ratio", "im_part"), file = data_filename)
  cat(">>> 模拟完成，数据已保存至", data_filename, "\n\n")
}

# 统计量计算模块
cat(">>> [2/4] 计算理论密度与统计量...\n")

# 定义斯蒂尔杰斯变换求解函数
solve_stieltjes <- function(z, a2, H) {
  K <- ncol(H); n <- nrow(H)
  beta_old <- rep(0 + 0i, K)
  for(i in 1:200) {
    denom <- apply(H, 1, function(lambda_i) sum(a2 * lambda_i * beta_old) + z)
    beta_new <- -sapply(1:K, function(r) mean(H[, r] / denom))
    if (mean(abs(beta_new - beta_old)) < 1e-8) break
    beta_old <- beta_new
  }
  denom <- apply(H, 1, function(lambda_i) sum(a2 * lambda_i * beta_old) + z)
  return(list(s = -mean(1 / denom)))
}

# 设定用于计算和绘图的 x 轴序列
x_seq_min <- min(min(eig_Cn), plot_xlim[1]) * 1.1
x_seq_max <- max(max(eig_Cn), plot_xlim[2]) * 1.1
x_seq <- seq(x_seq_min, x_seq_max, length.out = 300)

# 计算理论密度 (Density)
theo_density <- sapply(x_seq, function(x) {
  z <- complex(real = x, imaginary = im_part)
  return(Im(solve_stieltjes(z, a2, H)$s) / pi)
})

# 计算统计量 (KS距离 和 L2距离)
ecdf_fun <- ecdf(eig_Cn)
ecdf_vals <- ecdf_fun(x_seq)

dx <- diff(x_seq)[1]
cdf_theo_vals <- cumsum(theo_density) * dx
# 归一化 CDF 以消除数值积分误差的影响
if(max(cdf_theo_vals) > 0) cdf_theo_vals <- cdf_theo_vals / max(cdf_theo_vals)

ks_stat <- max(abs(ecdf_vals - cdf_theo_vals))
l2_stat <- sqrt(sum((ecdf_vals - cdf_theo_vals)^2) * dx)

cat("   -> KS Distance:", round(ks_stat, 4), "\n")
cat("   -> L2 Distance:", round(l2_stat, 4), "\n\n")


# 绘图
library(ggplot2)


if(!require(showtext)) install.packages("showtext")
library(showtext)
showtext_auto() 
font_add("SimSun", "simsun.ttc") 

if (!exists("K")) stop("变量 K, p, n 未定义，请先运行前面的模拟代码！")
plot_file_png <- paste0("plot_K", K, "_p", p, "_n", n, ".png")
plot_file_eps <- paste0("plot_K", K, "_p", p, "_n", n, ".eps")


# 直方图
hist_data <- data.frame(eigenvalue = eig_Cn)

# 理论曲线数据
theoretical_df <- data.frame(x = x_seq, density = theo_density)

custom_breaks <- seq(from = -3, to = 3, by = 0.05)

p_plot <- ggplot() +

  geom_histogram(data = hist_data, aes(x = eigenvalue, y = after_stat(density)),
                 breaks = custom_breaks,
                 fill = "#99CCFF",      
                 color = "white",       
                 alpha = 0.8,           
                 linewidth = 0.3) +     
  # 绘制理论曲线
  geom_line(data = theoretical_df, aes(x = x, y = density), 
            color = "#FF0000",           # 红色
            linewidth = 1.0) +           # 线宽
  
  coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
  
  labs(x = "特征值", y = "密度", 
       title = paste0("p=", p, ", n=", n, ", K=", K)) +
  
  theme_bw() + 
  theme(
    text = element_text(family = "SimSun"), 

    plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),

    axis.title = element_text(size = 26, face = "bold"),

    axis.text = element_text(size = 26, color = "black"),
    
    panel.grid = element_blank()
  )

#图片输出


cat(">>> [4/4] 生成 PNG 高清图...\n")
ggsave(plot_file_png, plot = p_plot, width = 8, height = 6, dpi = 300)
cat("PNG 已保存:", plot_file_png, "\n")

cat(">>> [额外] 生成 EPS 矢量图 (使用 Cairo 引擎)...\n")

ggsave(plot_file_eps, plot = p_plot, 
       device = cairo_ps, 
       width = 8, 
       height = 6, 
       fallback_resolution = 600)

cat("EPS 已保存:", plot_file_eps, "\n")
cat("提示：EPS 字体已锁定为 26pt，宋体，支持透明度。\n")

showtext_auto(FALSE)