library(statmod)

data=read.csv('C:\\Users\\mengxing\\Desktop\\审稿意见\\3\\敏感性分析\\expression_matrix.csv')
rownames(data) <- data[, 1]
data=data[, -1]

par_e=read.csv('C:\\Users\\mengxing\\Desktop\\审稿意见\\3\\敏感性分析\\par.csv')
rownames(par_e)=par_e[, 1]
par_e=par_e[, -1]


col_names=c('kon_moment1','kon_moment2','kon_moment3','kon_fano',
               'koff_moment1','koff_moment2','koff_moment3','koff_fano',
               'kb_moment1','kb_moment2','kb_moment3','kb_fano',
              'theta_moment1','theta_moment2','theta_moment3','theta_fano')
result = data.frame(matrix(NA, nrow=nrow(par_e), ncol=length(col_names)))
colnames(result) = col_names
row.names(result)=row.names(par_e)

BP_prob=function (n=NULL,kon=NULL, koff=NULL,ksyn=NULL,theta=NULL){

  fn3=function(x1, x2, m) {
    if (max(m) < 1e+05)
      #
      res=ppois(x2,m)-ppois(x1,m)
    #
    else res=pnorm(x2, m, sqrt(m)) - pnorm(x1, m, sqrt(m))
    return(res)
  }
  #
  w=gauss.quad(10, "jacobi", alpha=koff - 1, beta=kon - 1)
  gs=sum(w$weight * fn3(x1=n-1, x2=n, m=ksyn * (1 + w$node)/2))
  #
  if(n==0){
    prob=theta+1/beta(kon, koff) * 2^(-kon - koff + 1) * gs*(1-theta)
  }else{
    prob=1/beta(kon, koff) * 2^(-kon - koff + 1) * gs*(1 - theta)
  }
  return(prob)
}

ZI_telegraph_dist = function(kon, koff, kb, theta, N){
  pm = numeric(N+1)

  for(n in 0:N){
    pm[n+1] = BP_prob(n=n, kon=kon, koff=koff, ksyn=kb, theta=theta)
  }
  return(pm)
}

sfun_m1 = function(kon, koff, kb, theta, N){
  pm = ZI_telegraph_dist(kon, koff, kb, theta, N)
  sum( (0:N) * pm )
}

sfun_m2 = function(kon, koff, kb, theta, N){
  pm = ZI_telegraph_dist(kon, koff, kb, theta, N)
  sum( (0:N)^2 * pm )
}

sfun_m3 = function(kon, koff, kb, theta, N){
  pm = ZI_telegraph_dist(kon, koff, kb, theta, N)
  sum( (0:N)^3 * pm )
}

sfun_fano = function(kon, koff, kb, theta, N){
  pm = ZI_telegraph_dist(kon, koff, kb, theta, N)
  m1 = sum( (0:N) * pm )
  m2 = sum( (0:N)^2 * pm )
  (m2/m1) - m1
}

dt = 0.0000001

for(zushu in row.names(data)){
  pmdata = data[zushu,]
  pmdata = pmdata[!is.na(pmdata)]
  N=max(pmdata)

  kon_e  = par_e[zushu,]$kon_ctrl
  koff_e = par_e[zushu,]$koff_ctrl
  kb_e   = par_e[zushu,]$ksyn_ctrl
  theta_e = par_e[zushu,]$theta_ctrl

  s_kon_m1 = abs((sfun_m1(kon_e+dt,koff_e,kb_e,theta_e,N) - sfun_m1(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kon_m2 = abs((sfun_m2(kon_e+dt,koff_e,kb_e,theta_e,N) - sfun_m2(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kon_m3 = abs((sfun_m3(kon_e+dt,koff_e,kb_e,theta_e,N) - sfun_m3(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kon_fano = abs((sfun_fano(kon_e+dt,koff_e,kb_e,theta_e,N)-sfun_fano(kon_e,koff_e,kb_e,theta_e,N))/dt)

  s_koff_m1 = abs((sfun_m1(kon_e,koff_e+dt,kb_e,theta_e,N) - sfun_m1(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_koff_m2 = abs((sfun_m2(kon_e,koff_e+dt,kb_e,theta_e,N) - sfun_m2(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_koff_m3 = abs((sfun_m3(kon_e,koff_e+dt,kb_e,theta_e,N) - sfun_m3(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_koff_fano = abs((sfun_fano(kon_e,koff_e+dt,kb_e,theta_e,N)-sfun_fano(kon_e,koff_e,kb_e,theta_e,N))/dt)

  s_kb_m1 = abs((sfun_m1(kon_e,koff_e,kb_e+dt,theta_e,N) - sfun_m1(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kb_m2 = abs((sfun_m2(kon_e,koff_e,kb_e+dt,theta_e,N) - sfun_m2(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kb_m3 = abs((sfun_m3(kon_e,koff_e,kb_e+dt,theta_e,N) - sfun_m3(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_kb_fano = abs((sfun_fano(kon_e,koff_e,kb_e+dt,theta_e,N)-sfun_fano(kon_e,koff_e,kb_e,theta_e,N))/dt)

  s_theta_m1 = abs((sfun_m1(kon_e,koff_e,kb_e,theta_e+dt,N) - sfun_m1(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_theta_m2 = abs((sfun_m2(kon_e,koff_e,kb_e,theta_e+dt,N) - sfun_m2(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_theta_m3 = abs((sfun_m3(kon_e,koff_e,kb_e,theta_e+dt,N) - sfun_m3(kon_e,koff_e,kb_e,theta_e,N))/dt)
  s_theta_fano = abs((sfun_fano(kon_e,koff_e,kb_e,theta_e+dt,N)-sfun_fano(kon_e,koff_e,kb_e,theta_e,N))/dt)

  result[zushu,] = c(
    s_kon_m1, s_kon_m2, s_kon_m3, s_kon_fano,
    s_koff_m1, s_koff_m2, s_koff_m3, s_koff_fano,
    s_kb_m1, s_kb_m2, s_kb_m3, s_kb_fano,
    s_theta_m1, s_theta_m2, s_theta_m3, s_theta_fano
  )

  cat("完成基因：", zushu, "\n")
}

library(tidyr)
library(dplyr)
library(ggplot2)
write.csv(x=result, file='result_FSP_Sensitivity_moment.csv',row.names=TRUE)

result$Gene <- rownames(result)

df_long <- result %>%
  pivot_longer(
    cols = -Gene,
    names_to = c("Parameter","Statistic"),
    names_sep = "_",
    values_to = "Sensitivity"
  ) %>%
  mutate(
    Parameter = factor(Parameter, levels = c("kon","koff","kb","theta")),

    Statistic = factor(Statistic, levels = c("moment1","moment2","moment3"))
  )


df_clean <- df_long %>%
  group_by(Parameter, Statistic) %>%
  filter(Sensitivity < quantile(Sensitivity, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Parameter, Statistic) %>%
  filter(Sensitivity < quantile(Sensitivity, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(log10_Sensitivity = log10(Sensitivity))


ggplot(df_clean, aes(x = Statistic, y = log10_Sensitivity, color = Parameter)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  labs(x = "Statistical index", y = "log10(Sensitivity)", color = "Parameters") +
  theme_bw() +
  ylim(-2, 3)



