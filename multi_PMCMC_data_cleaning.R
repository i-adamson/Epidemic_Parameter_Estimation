s <- multi_bi_lst$S[,3]
i1 <- multi_bi_lst$I1[,3]
i2 <- multi_bi_lst$I2[,3]
g <- multi_bi_lst$gamma
a <- multi_bi_lst$alpha

time_deps <- bind_cols(beta_1, beta_2[,-(1:2)],s,i1,i2 )
colnames(time_deps)<- c('np', 'time', 'beta1', 'beta2', 'S', 'I1', 'I2')


beta_df <- betas %>%
  group_by(np)%>%
  group_split()
View(beta_df)


aux <- bind_cols(g,a[-1])
colnames(aux) <- c('np', 'gamma','alpha')


write.csv(time_deps, "C:/Users/nwwt55/Documents/post_dat.csv", row.names = FALSE)
write.csv(aux, "C:/Users/nwwt55/Documents/aux_params.csv", row.names = FALSE)
