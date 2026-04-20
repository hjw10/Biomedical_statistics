setwd("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework3/")

## question1
# (1)
# load data
data(mtcars)
mtcars
# grep data
data_4cyl = mtcars[mtcars$cyl == 4,]

# check the necessary assumptions
shapiro.test(data_4cyl$mpg)
# normal 

# t test
t.test(data_4cyl$mpg, mu = 23, alternative = "greater")
# the mean miles per gallon (mpg) of cars with 4 cylinders is significantly greater than 23 mpg.

# (2)
set.seed(123)  # 为了结果可重复
observed_mean <- mean(data_4cyl$mpg)
shifted_sample <- data_4cyl$mpg - observed_mean + 23

n_perm <- 10000  # 置换次数
boot_means <- replicate(n_perm, mean(sample(shifted_sample, replace = TRUE)))

# 右侧单侧检验：计算重抽样均值 >= 观测均值的比例
p_perm <- mean(boot_means >= observed_mean)
p_perm # 0.0028


# question2
score1 = c(75,80,70,85,65,90,78,72,88,76,84,70,68,85,77)
score2 = c(85,82,78,88,70,95,84,80,92,81,87,73,75,90,83)
score = as.data.frame(cbind(score1,score2))
colnames(score) = c("Before", "After")
# t test
t.test(score$After, score$Before, paired = T, conf.level = 0.99, alternative = 'greater')

# fig
boxplot(score$Before, score$After,
        names = c("Before", "After"),
        ylab = "score",
        col = c("lightblue", "lightcoral"))


# question3
Cr_placebo   = c(0.00, 0.01, -0.03, -0.01, -0.02, -0.01, 0.03, -0.02, 0.00, 0.04, 0.02, -0.02, -0.01, 0.00, -0.03)
Cr_treatment = c(0.23, 0.08,  0.10,  0.15,  0.09,  0.12, 0.19,  0.08, 0.16, 0.10, 0.19,  0.12,  0.15, 0.13,  0.14)

var.test(Cr_treatment, Cr_placebo)
t.test(Cr_treatment, Cr_placebo, alternative = "greater", var.equal = F)


# question4
independent_t_test <- function(group1,group2,alpha = 0.05) {
  mean1 = mean(group1)
  mean2 = mean(group2)
  mean_diff = mean1 - mean2
  n1 = length(group1)
  n2 = length(group2)
  var1 = var(group1)
  var2 = var(group2)
  pooled_var <- ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
  se <- sqrt(pooled_var * (1/n1 + 1/n2))
  t_stat <- mean_diff / se
  df <- n1 + n2 - 2
  t_crit <- qt(1 - alpha/2, df)
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  if (abs(t_stat) > t_crit) {
    conclusion <- paste("拒绝 H0：两组均值存在显著差异 (p <", alpha, ")")
  } else {
    conclusion <- paste("不拒绝 H0：两组均值无显著差异 (p >=", alpha, ")")
  }
  result <- list(
    sample1 = list( mean = mean1, variance = var1),
    sample2 = list( mean = mean2, variance = var2),
    mean_difference = mean_diff,
    pooled_variance = pooled_var,
    standard_error = se,
    t_statistic = t_stat,
    df = df,
    critical_value = t_crit,
    p_value = p_value,
    conclusion = conclusion
  )
  return(result)
}

group_a = c(18.5, 19.2, 17.8, 20.1, 19.5, 18.9)
group_b = c(15.2, 16.8, 14.9, 16.3, 15.8, 16.1)

result = independent_t_test(group_a, group_b)
result


# question5
# (1) 配对检验；数据非正态分布、小样本
# (2)
viral_data <- data.frame(
  Sample_ID = 1:10,
  Before = c(6.5, 7.2, 5.8, 8.1, 6.9, 7.5, 6.3, 7.8, 7.0, 6.7),
  After = c(5.8, 7.0, 5.6, 7.2, 6.1, 7.6, 5.5, 7.1, 6.8, 6.2)
)
wilcox_result <- wilcox.test(viral_data$After, viral_data$Before,
                             paired = TRUE,
                             alternative = "two.sided",
                             exact = FALSE,
                             correct = TRUE)
wilcox_result
# (3)
boxplot(viral_data$Before, viral_data$After,
        names = c("Before", "After"),
        ylab = "log10 viral loads",
        col = c("lightblue", "lightcoral"))

# (4)零假设认为治疗前后样本的病毒载量不变，
# p-value小于0.01表示零假设成立的条件下，我们得到现在的实验数据的概率小于0.01，
# 这种可能性很小，因此我们拒绝零假设，采纳备择假设认为治疗前后的病毒载量改变




