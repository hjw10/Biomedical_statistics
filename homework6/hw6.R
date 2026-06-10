setwd("C:/Users/hujw0/Desktop/学习资料/生物医学统计/homework6/")

# Q1
data(PlantGrowth)
str(PlantGrowth)
head(PlantGrowth)

# 拆分数据
trt1_weights <- PlantGrowth$weight[PlantGrowth$group == "trt1"]
trt2_weights <- PlantGrowth$weight[PlantGrowth$group == "trt2"]
pg_sub <- PlantGrowth[PlantGrowth$group %in% c("trt1", "trt2"), ]
# Shapiro-Wilk 检验
shapiro.test(trt1_weights)
shapiro.test(trt2_weights)

# 方差齐性检验
library(car)
leveneTest(weight ~ group, data = pg_sub)

# 显著性
t_test_result <- t.test(weight ~ group, data = pg_sub, var.equal = TRUE)
t_test_result

# Q2
data2 = read.csv("homework6-Q2.csv", row.names = 1)
# 单因素 ANOVA
aov1 <- aov(viability ~ drug, data = data2)
summary(aov1)
# 进行 Tukey HSD 事后检验
TukeyHSD(aov1)

# 双因素 ANOVA
aov2 <- aov(viability ~ drug * time, data = data2)
summary(aov2)
# 创建组合因子
data2$comb <- interaction(data2$drug, data2$time)
# 两两 t 检验
pairwise.t.test(data2$viability, data2$comb, p.adjust.method = "BH", pool.sd = FALSE)


# Q3
# 1. 散点图矩阵和线性回归
data(swiss)
df <- swiss[, c("Fertility", "Education", "Agriculture", "Infant.Mortality")]
pairs(df)
model_full <- lm(Fertility ~ Education + Agriculture + Infant.Mortality, data = df)
summary(model_full)

# 2. 简化模型
summary(model_full)$r.squared
summary(model_full)$adj.r.squared
# 简化模型去掉 Agriculture
model_reduced <- lm(Fertility ~ Education + Infant.Mortality, data = df)
anova(model_reduced, model_full) 

# 3. 逻辑回归
df$high_fertility <- ifelse(df$Fertility > median(df$Fertility), 1, 0)
logit_model <- glm(high_fertility ~ Education + Agriculture + Infant.Mortality, 
                   data = df, family = binomial)
summary(logit_model)

# Q4
data(USArrests)
pca <- prcomp(USArrests, scale. = TRUE)   # 标准化后 PCA
summary(pca) # 3

biplot(pca, main = "PCA biplot of USArrests")

scores <- pca$x[, 1:2]          # 提取前两个主成分得分
set.seed(123)                    # 固定随机种子
km <- kmeans(scores, centers = 3, nstart = 25)
cluster_assign <- km$cluster
names(cluster_assign) <- rownames(USArrests)
print(cluster_assign)
aggregate(USArrests, by = list(Cluster = km$cluster), FUN = mean)

library(ggplot2)
plot_df <- data.frame(PC1 = scores[,1], PC2 = scores[,2], Cluster = factor(km$cluster))
ggplot(plot_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering on PC1-PC2 Space (USArrests)") +
  theme_minimal()

# Q5





