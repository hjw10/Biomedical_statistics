setwd("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework4/")
library(car)

## question1
# (1)
data = read.csv("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework4/homework4-dataset/homework4-Q1-data.csv", row.names = 1)
data <- as.data.frame(t(data))  # 转置：样本为行，基因为列
data[] <- lapply(data, as.numeric)

# 定义分组
group <- ifelse(grepl("^S", rownames(data)), "Sensitive", "Resistant")
S_idx <- which(group == "Sensitive")
R_idx <- which(group == "Resistant")

# 存储结果
res <- data.frame(
  Gene      = colnames(data),
  Method    = character(12),
  p_value   = numeric(12),
  log2FC    = numeric(12),
  stringsAsFactors = FALSE
)

for(i in 1:ncol(data)) {
  gene <- colnames(data)[i]
  s <- data[S_idx, i]
  r <- data[R_idx, i]
  
  # Shapiro-Wilk 正态性检验
  p.norm.S <- shapiro.test(s)$p.value
  p.norm.R <- shapiro.test(r)$p.value
  normal <- (p.norm.S > 0.05) && (p.norm.R > 0.05)
  
  equal_var <- FALSE
  if(normal) {
    y <- c(s, r)
    grp <- factor(c(rep("S",10), rep("R",10)))
    p.var <- leveneTest(y ~ grp)$`Pr(>F)`[1]
    equal_var <- (p.var > 0.05)
  }
  
  # 选择检验方法
  if (normal && equal_var) {
    meth <- "Student t-test"
    pval <- t.test(s, r, var.equal = TRUE)$p.value
  } else if (normal && !equal_var) {
    meth <- "Welch t-test"
    pval <- t.test(s, r, var.equal = FALSE)$p.value
  } else {
    meth <- "Wilcoxon rank-sum test"
    pval <- wilcox.test(s, r, exact = FALSE)$p.value
  }
  
  log2FC <- log2(mean(r) / mean(s))
  
  res[i, "Method"] <- meth
  res[i, "p_value"] <- pval
  res[i, "log2FC"] <- log2FC
}

# BH 校正
res$p_adj_BH <- p.adjust(res$p_value, method = "BH")
print(res)



# Question2
library(pwr)

# (1) 给定每组 n = 35，d = 0.6，α = 0.05，求 power
power1 <- pwr.t.test(n = 35, d = 0.6, sig.level = 0.05,
                     type = "two.sample", alternative = "two.sided")
power1$power

# (2) 给定 d = 0.6，α = 0.05，power = 0.8，求每组所需样本量 n
n_needed <- pwr.t.test(d = 0.6, sig.level = 0.05, power = 0.8,
                       type = "two.sample", alternative = "two.sided")
n_needed$n

# (3) 给定每组 n = 20，α = 0.05，power = 0.8，求最小可检测效应量 d
d_detect <- pwr.t.test(n = 20, power = 0.8, sig.level = 0.05,
                       type = "two.sample", alternative = "two.sided")
d_detect$d



# Question3
# (1)
data = read.csv("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework4/homework4-dataset/homework4-Q3-data.csv")
data$CellLine <- as.factor(data$CellLine)

boxplot(GFP_Positive_Percent ~ CellLine, data = data,
        main = "Infection Efficiency by Cell Line",
        xlab = "Cell Line", ylab = "% GFP-positive cells",
        col = c("gold", "tomato", "skyblue", "lightgreen"),
        las = 1)

tapply(data$GFP_Positive_Percent, data$CellLine, mean)
tapply(data$GFP_Positive_Percent, data$CellLine, sd)

anova_result <- aov(GFP_Positive_Percent ~ CellLine, data = data)
summary(anova_result)            # 展示 ANOVA 表

# 从 summary() 中提取 p 值
p_anova <- summary(anova_result)[[1]]$`Pr(>F)`[1]

# (2)
if (p_anova < 0.05) {
  cat("ANOVA is significant. Performing Tukey's HSD...\n")
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  # 可视化 HSD 结果
  plot(tukey_result, las = 1)
} else {
  cat("ANOVA not significant. No post-hoc test performed.\n")
}


# Question4
data = read.csv("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework4/homework4-dataset/homework4-Q4-data.csv")
data$A <- factor(data$A)
data$B <- factor(data$B)
data$group <- interaction(data$A, data$B)
# (1)
model_full <- aov(Weight ~ A * B, data = data)
shapiro.test(residuals(model_full))

leveneTest(Weight ~ A * B, data = data)
# (2)
# 主效应模型
model_main <- aov(Weight ~ A + B, data = data)
summary(model_main)
# (3)
model_inter <- aov(Weight ~ A * B, data = data)
summary(model_inter)

library(emmeans)
# 基于交互效应模型
emm <- emmeans(model_inter, ~ A * B)
# 以 A2B1 为参考组，进行 Dunnett 校正
contrast(emm, method = "dunnett", ref = "A2 B1")


# Question5
drugA <- c(145, 152, 138, 160, 147, 155, 142, 149, 153, 140)
drugB <- c(136, 129, 142, 131, 140, 135, 138, 133, 145, 130)
drugC <- c(120, 125, 118, 132, 121, 128, 115, 124, 130, 122)

weight <- c(drugA, drugB, drugC)
group <- factor(rep(c("A", "B", "C"), each = 10))
data <- data.frame(weight, group)

# 由于数据不满足正态性和方差齐性，使用 Kruskal-Wallis 检验
kruskal.test(weight ~ group, data = data)

library(FSA)
# Dunn 检验，双侧，Bonferroni 校正
dunnTest(weight ~ group, data = data, method = "bonferroni")



