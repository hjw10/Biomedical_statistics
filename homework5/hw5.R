setwd("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework5/")
library(ggplot2)


# Question 1
data1 = read.table("homework5-Q1.tsv", sep = '\t', header = TRUE)
rownames(data1) = data1$GeneID
data1 = data1[,-1]

# (1)
GAPDH_expr <- as.numeric(data1["GAPDH", ])
ACTB_expr  <- as.numeric(data1["ACTB", ])


# pearson
pearson_manual <- function(x, y) {
  n <- length(x)
  # 去除缺失值（若有）
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  n <- length(x)
  mean_x <- mean(x)
  mean_y <- mean(y)
  # 计算分子（协方差）和分母（标准差乘积）
  num <- sum((x - mean_x) * (y - mean_y))
  denom <- sqrt(sum((x - mean_x)^2) * sum((y - mean_y)^2))
  r <- num / denom
  return(r)
}

pearson_r <- pearson_manual(GAPDH_expr, ACTB_expr)
print(paste("Pearson r =", round(pearson_r, 4)))


# spearman
spearman_manual <- function(x, y) {
  n <- length(x)
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  n <- length(x)
  # 计算秩
  rank_x <- rank(x, ties.method = "average")
  rank_y <- rank(y, ties.method = "average")
  # 对秩计算 Pearson 相关系数
  mean_rank_x <- mean(rank_x)
  mean_rank_y <- mean(rank_y)
  num <- sum((rank_x - mean_rank_x) * (rank_y - mean_rank_y))
  denom <- sqrt(sum((rank_x - mean_rank_x)^2) * sum((rank_y - mean_rank_y)^2))
  
  rho <- num / denom
  return(rho)
}

spearman_rho <- spearman_manual(GAPDH_expr, ACTB_expr)
print(paste("Spearman rho =", round(spearman_rho, 4)))


# Question 2
age <- c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 33, 43, 53, 63, 73)
sbp <- c(115, 122, 128, 132, 140, 148, 155, 160, 168, 175, 118, 124, 130, 135, 142, 150, 157, 162, 170, 178, 120, 131, 140, 151, 161)
bp <- data.frame(age, sbp)

# 线性拟合
model <- lm(sbp ~ age, data = bp)
summary(model)
# 回归函数
coef(model)
# 95% 置信区间
conf_int <- confint(model, "age", level = 0.95)
conf_int


# Question 3
data(mtcars)
# （1）构建模型
model <- lm(mpg ~ wt + hp, data = mtcars)
summary(model)
# （2）提取p value
coef_summary <- summary(model)$coefficients
p_wt <- coef_summary["wt", "Pr(>|t|)"]
p_hp <- coef_summary["hp", "Pr(>|t|)"]

cat("wt 系数的 p 值 =", p_wt, "\n")
cat("hp 系数的 p 值 =", p_hp, "\n")

# （3） 提取F检验p value
f_p_value <- summary(model)$fstatistic
p_value_f <- pf(f_p_value[1], f_p_value[2], f_p_value[3], lower.tail = FALSE)


# Question 4
Age <- c(1,2,3,4,5,6,7,8)
Hormone_Ratio <- c(0.42, 0.39, 0.53, 0.65, 0.75, 0.56, 0.53, 0.48)
# （1）
# 计算自然对数
ln_ratio <- log(Hormone_Ratio)
# 散点图
plot(Age, ln_ratio, 
     xlab = "Age (years)", ylab = "ln(Hormone Ratio)",
     main = "ln(Hormone Ratio) vs Age", pch = 19, col = "blue")

# （2）
model_linear <- lm(ln_ratio ~ Age)
summary(model_linear)
# 在原图上添加回归线
abline(model_linear, col = "red", lwd = 2)

# （3）
# 二次模型
Age2 <- Age^2
model_quad <- lm(ln_ratio ~ Age + Age2)
summary(model_quad)
# 绘制二次曲线（预测值连线）
Age_seq <- seq(min(Age), max(Age), length.out = 100)
pred_quad <- predict(model_quad, newdata = data.frame(Age = Age_seq, Age2 = Age_seq^2))
lines(Age_seq, pred_quad, col = "green", lwd = 2)
# 方差分析：比较线性模型与二次模型
anova(model_linear, model_quad)


# Question 5
library(MASS)
data(birthwt)
birthwt$race <- factor(birthwt$race, labels = c("white", "black", "other"))
birthwt$smoke <- factor(birthwt$smoke, labels = c("non-smoker", "smoker"))
birthwt$ht <- factor(birthwt$ht, labels = c("no", "yes"))
str(birthwt)

# a
# 完整模型
model_full <- glm(low ~ age + lwt + race + smoke + ht, 
                  data = birthwt, family = binomial)
summary(model_full)

# b
# 检验 race 整体的显著性
anova(model_full, test = "Chisq")
# 简化模型
model_reduced <- glm(low ~ lwt + smoke + ht, data = birthwt, family = binomial)
summary(model_reduced)

# c
anova(model_reduced, model_full, test = "Chisq")

# d
set.seed(320)

# 数据集共 189 行
n <- nrow(birthwt)
train_idx <- sample(1:n, size = 100, replace = FALSE)
train_data <- birthwt[train_idx, ]
test_data <- birthwt[-train_idx, ]

# 用所有变量 (age, lwt, race, smoke, ht) 拟合训练集模型
model_train <- glm(low ~ age + lwt + race + smoke + ht, 
                   data = train_data, family = binomial)

# 预测测试集概率
pred_prob <- predict(model_train, newdata = test_data, type = "response")
pred_class <- ifelse(pred_prob >= 0.5, 1, 0)

# 真实类别
true_class <- test_data$low

accuracy <- mean(pred_class == true_class)
cat("预测准确率 =", round(accuracy, 4), "\n")

# 也可输出混淆矩阵
table(pred_class, true_class)









