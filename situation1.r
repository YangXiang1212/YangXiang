
rm(list = ls()) # 清空环境
library(ggplot2)

if(!require(showtext)) install.packages("showtext")
library(showtext)

showtext_auto()
font_add("SimSun", "simsun.ttc") 

#函数定义
simulate_sn_eigvals <- function(params) {
  K <- params$K
  p <- params$p
  n <- params$n
  a <- params$a
  b <- params$b
  
  X_n <- matrix(complex(real = 0, imaginary = 0), nrow = p, ncol = n)
  Tr <- 0  
  
  for (r in 1:K) {
    Z_real <- matrix(rnorm(p * n), p, n)
    Z_imag <- matrix(rnorm(p * n), p, n)
    Z_r <- (Z_real + 1i * Z_imag) / sqrt(2)  
    
    X_n <- X_n + sqrt(a[r]) * Z_r * sqrt(b[r])  
    Tr <- Tr + a[r] * b[r]  
  }
  
  term1 <- (Conj(t(X_n)) %*% X_n) / p  
  term2 <- Tr * diag(n)  
  S <- (term1 + Conj(t(term1))) / 2  
  S_n <- sqrt(p / n) * (S - term2)  
  
  S_n_clean <- zapsmall(Re(S_n)) 
  S_n_clean[is.infinite(S_n_clean) | is.na(S_n_clean)] <- 0 
  eigvals <- eigen(S_n_clean, symmetric = TRUE, only.values = TRUE)$values
  return(eigvals)
}

stieltjes_s <- function(z, gamma) {
  sqrt_term <- sqrt(z^2 - 4 * gamma)
  if (Re(z) < 0 && Im(sqrt_term) < 0) {
    sqrt_term <- -sqrt_term
  }
  s <- -1/(2 * gamma) * (z - sqrt_term)
  return(s)
}

inverse_stieltjes <- function(x, gamma, eps = 1e-4) {
  z <- x + eps * 1i
  s <- stieltjes_s(z, gamma)
  density <- Im(s) / pi
  return(max(density, 0))
}

#参数设置
K <- 2 
a <- runif(K, min = 1, max = 3) 
b <- runif(K, min = 1, max = 3) 

params <- list(
  K = K,
  p = 150000,    
  n = 1000,    
  a = a,
  b = b
)

file_prefix <- paste0("K", K, "p", params$p, "n", params$n)
data_filename <- paste0(file_prefix, ".RData")

#数据检测与加载
if (file.exists(data_filename)) {
  cat("✅ 检测到已保存数据文件，直接加载...\n")
  load(data_filename)
  gamma <- save_data$gamma
  eigvals <- save_data$eigvals
  theoretical_density <- save_data$theoretical_density
  # 重新定义 x_range 以防万一
  theo_left <- -2 * sqrt(gamma)
  theo_right <- 2 * sqrt(gamma)
  x_range <- save_data$x_range
  
} else {
  cat("❌ 未检测到数据文件，开始执行模拟...\n")



gamma <- sum(params$a * (params$b)^2)

cat("Adjusted gamma for Complex case =", round(gamma, 2), "\n")

  eigvals <- simulate_sn_eigvals(params)
  
  theo_left <- -2 * sqrt(gamma)
  theo_right <- 2 * sqrt(gamma)
  
  x_min <- theo_left - 0.5
  x_max <- theo_right + 0.5
  x_range <- seq(x_min, x_max, by = 0.1)
  
  theoretical_density <- sapply(x_range, inverse_stieltjes, gamma = gamma)
  
  # 保存数据
  save_data <- list(params = params, gamma = gamma, eigvals = eigvals,
                    theoretical_density = theoretical_density, x_range = x_range)
  save(save_data, file = data_filename, compress = TRUE)
}

#绘图

theoretical_df <- data.frame(x = x_range, density = theoretical_density)
hist_data <- data.frame(eigenvalue = eigvals)

plot <- ggplot() +

  geom_histogram(data = hist_data, aes(x = eigenvalue, y = after_stat(density)),
                 bins = 50, 
                 fill = "#99CCFF", 
                 color = "white",  
                 alpha = 0.8,      
                 linewidth = 0.3) +
  
  # 理论曲线
  geom_line(data = theoretical_df, aes(x = x, y = density), 
            color = "#FF0000", linewidth = 1.0) +
  
  xlim(theo_left - 0.5, theo_right + 0.5) +
  
  labs(x = "特征值", y = "密度", 
       title = paste0("K=", K, ", gamma=", round(gamma, 2))) +
  
  theme_bw() +
  #中文字体
  theme(
    text = element_text(family = "SimSun"),
    plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
    axis.title = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 26),
    panel.grid = element_blank() 
  )

# 图片保存

png_filename <- paste0(file_prefix, ".png")
ggsave(png_filename, plot, width = 8, height = 6, dpi = 300)
cat("PNG 已保存:", png_filename, "\n")


eps_filename <- paste0(file_prefix, ".eps")

ggsave(eps_filename, plot, device = cairo_ps, width = 8, height = 6, fallback_resolution = 600)

cat("EPS 已保存:", eps_filename, "\n")
cat("提示：此 EPS 文件使用 showtext 生成，文字已转为矢量路径，放入 LaTeX 绝不会乱码。\n")


showtext_auto(FALSE)