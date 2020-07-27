library(knockoff)
library("doParallel")      
library("foreach")  


#model_parameters
n = 1000
q = 0.2

nonzero = sample(d,sp)
beta = numeric(d)

agg = function(k,set_list,c){
    m = numeric(d)
    ss = numeric(k)
    for(i in 1:k){
        s = set_list[[i]]
        ss[i] = length(s)
        for(j in 1:d){
            if(sum(s == j)==1){
                m[j] = m[j] + 1
            }
        }
    }
    set_agg = (1:d)[m >= c]
    return(list(set = set_agg, m = m))
}

p_values = function(stat){
    d = length(stat)
    p_v = numeric(d)
    for(j in 1:d){
        if(stat[j] <= 0){
            p_v[j] = 1
        }
        else{
            p_v[j] = (1/d)*(1+sum(stat <= -stat[j]))
        }
    }
    return(p_v)
}

AKO = function(pval_list,d,k){
    gamma = 0.5
    p_agg = numeric(d)
    for(j in 1:d){
        pval_b = numeric(k)
        for(i in 1:k){
            pval_b[i] = pval_list[[i]][j]
        }
        p_agg = min(1,(1/gamma)*quantile(pval_b,gamma))
    }
    p_agg_rank = p_agg[order(p_agg)]
    kagg = max((p_agg_rank <= (log(d))*(1:d)*q/d) == T)
    return(AKO_select = (1:d)[p_agg <= p_agg_rank[kagg]])
}

agg_adp = function(k,set_list){
    if(k==1) return(list(set = set_list[[1]],
                         c = 1,cset = 1,
                         ratio = 1,size = length(set_list)))
    m = numeric(d)
    ss = numeric(k)
    for(i in 1:k){
        s = set_list[[i]]
        ss[i] = length(s)
        for(j in 1:d){
            if(sum(s == j)==1){
                m[j] = m[j] + 1
            }
        }
    }
    M = max(ss)
    mn = mean(ss)
    ssz = numeric(k)
    ratio = numeric(k-1)
    cc = NULL
    for(c in 1:k){
        set_agg = (1:d)[m >= c]
        ssz[c] = length(set_agg)
        if(c > 1) ratio[c-1] = 
            (ssz[c-1]+1)/(ssz[c]+1)
        if(length(set_agg) >= mn){
            cc = c(cc,c)
        }
    }
    ll = min(k-1,length(cc))
    #print(ratio[1:length(cc)])
    if(k <= 2) c_op = k
    if(k > 2) c_op = max(cc[ratio[1:ll] == min(ratio[1:ll])])
    set_agg = (1:d)[m >= c_op]
    return(list(set = set_agg,c = c_op,cset = cc,
                ratio = ratio,size = ssz,m = m))
}

agg_min = function(k,set_list){
    if(k==1) return(list(set = set_list[[1]],
                         c = 1,cset = 1,
                         ratio = 1,size = length(set_list)))
    m = numeric(d)
    ss = numeric(k)
    for(i in 1:k){
        s = set_list[[i]]
        ss[i] = length(s)
        for(j in 1:d){
            if(sum(s == j)==1){
                m[j] = m[j] + 1
            }
        }
    }
    M = max(ss)
    mn = mean(ss)
    ssz = numeric(k)
    ratio = numeric(k-1)
    cc = NULL
    for(c in 1:k){
        set_agg = (1:d)[m >= c]
        ssz[c] = length(set_agg)
    }
    obj = ssz * (1:k)
    c_op = (1:k)[obj == min(obj)]
    set_min = (1:d)[m >= c_op]
    return(list(set = set_min,c = c_op,
                size = ssz))
}

simu_i = function(d,k,sp){
    library(knockoff)
    nonzero = sample(d,sp)
    beta = numeric(d) 
    indv_set = NULL
    p_val = NULL
    indv_set1 = NULL
    nn = n/k
    beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
    for(i in 1:k){
        mu = rep(0,d)
        rho = 0.25
        Sigma = toeplitz(rho^(0:(d-1)))
        X = matrix(rnorm(nn*d),nn) %*% chol(Sigma)
        y.sample = function(X) X %*% beta + rnorm(nn)
        y = y.sample(X)
        result = knockoff.filter(X, y, 
                                 knockoffs = create.second_order, 
                                 statistic = stat.glmnet_lambdasmax,
                                 fdr = q)
        selected_set = result$selected
        indv_set[[i]] = result$selected
        # i_fdp = sum(beta[selected_set] == 0) / max(1, length(selected))
        # i_pw = sum(beta[selected_set] != 0) / length(nonzero)
        # ifdp[rn,i] = i_fdp
        # ipw[rn,i] = i_pw
        result1 = knockoff.filter(X, y, 
                                 knockoffs = create.second_order, 
                                 statistic = stat.glmnet_lambdasmax,
                                 fdr = q/k)
        selected_set1 = result1$selected
        indv_set1[[i]] = selected_set1
        p_val[[i]] = p_values(result$statistic)
    }
    
    #result_addaptive
    result_adp = agg_adp(k,indv_set)
    c = result_adp$c
    cset = result_adp$cset
    selected = result_adp$set
    FDP = sum(beta[selected] == 0) / max(1, length(selected))
    PW = sum(beta[selected] != 0) / length(nonzero)
    
    #result_min
    result_min = agg_min(k,indv_set)
    c_min = result_min$c
    selected = result_min$set
    FDP_min = sum(beta[selected] == 0) / max(1, length(selected))
    PW_min = sum(beta[selected] != 0) / length(nonzero)
    
    #AKO
    AKO_select = AKO(p_val, d, k)
    Afdp = sum(beta[AKO_select] == 0) / max(1, length(AKO_select))
    Apw = sum(beta[AKO_select] != 0) / length(nonzero)
    
    #union
    kk = 1
    result = agg(k,indv_set,kk)
    selected = result$set
    ufdp = sum(beta[selected] == 0) / max(1, length(selected))
    upw = sum(beta[selected] != 0) / length(nonzero)
    
    #intersection
    kk = k
    result = agg(k,indv_set,kk)
    selected = result$set
    ifdp = sum(beta[selected] == 0) / max(1, length(selected))
    ipw = sum(beta[selected] != 0) / length(nonzero)
    
    #Median
    kk = floor((k+1)/2)
    result = agg(k,indv_set,kk)
    selected = result$set
    mfdp = sum(beta[selected] == 0) / max(1, length(selected))
    mpw = sum(beta[selected] != 0) / length(nonzero)
    
    #Xie
    kk = 1
    result = agg(k,indv_set1,kk)
    selected = result$set
    Xfdp = sum(beta[selected] == 0) / max(1, length(selected))
    Xpw = sum(beta[selected] != 0) / length(nonzero)
    return(c(FDP,FDP_min,ifdp,ufdp,mfdp,Xfdp,Afdp,
                 PW,PW_min,ipw,upw,mpw,Xpw,Apw))
}




# FDP = numeric(repN)
# PW = numeric(repN)
# ifdp = matrix(rep(0,repN*k),nrow = repN)
# ipw = matrix(rep(0,repN*k),nrow = repN)
# fdp_seq = numeric(k)
# pw_seq = numeric(k)
# ssize = numeric(k)
# for(rn = 1:repN){
#     #distributed generation
#     indv_set = NULL
#     p_val = NULL
#     nn = n/k
#     beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
#     for(i in 1:k){
#         mu = rep(0,d)
#         rho = 0.25
#         Sigma = toeplitz(rho^(0:(d-1)))
#         X = matrix(rnorm(nn*d),nn) %*% chol(Sigma)
#         y.sample = function(X) X %*% beta + rnorm(nn)
#         y = y.sample(X)
#         result = knockoff.filter(X, y, 
#                                  knockoffs = create.fixed, 
#                                  statistic = stat.glmnet_lambdasmax,
#                                  fdr = q)
#         selected_set = result$selected
#         indv_set[[i]] = result$selected
#         p_val[[i]] = p_values(result$statistic)
#         i_fdp = sum(beta[selected_set] == 0) / max(1, length(selected))
#         i_pw = sum(beta[selected_set] != 0) / length(nonzero)
#         ifdp[rn,i] = i_fdp
#         ipw[rn,i] = i_pw
#     }
#     
#     #result_addaptive
#     result_adp = agg_adp(k,indv_set)
#     c = result_adp$c
#     cset = result_adp$cset
#     selected = result_adp$set
#     fdp = sum(beta[selected] == 0) / max(1, length(selected))
#     pw = sum(beta[selected] != 0) / length(nonzero)
#     #print(paste0("FDP = ",fdp,", Power = ",pw))
#     FDP[rn] = fdp
#     PW[rn] = pw
#     
#     #fixed_threshold
#     for(kk in 1:k){
#         result = agg(k,indv_set,kk)
#         selected = result$set
#         fdp = sum(beta[selected] == 0) / max(1, length(selected))
#         pw = sum(beta[selected] != 0) / length(nonzero)
#         fdp_seq[kk] = fdp + fdp_seq[kk]
#         pw_seq[kk] = pw + pw_seq[kk]
#         ssize[kk] = length(selected)
#     }
#     
# }
# 
# fdp_seq = fdp_seq/repN
# pw_seq = pw_seq/repN
# 
# fdpi = fdp_seq[k]
# pwi = pw_seq[k]
# fdpu = fdp_seq[1]
# pwu = pw_seq[1]
# 
# fdp_adp = mean(FDP)
# pw_adp = mean(PW)
# 
# par(mfrow=c(1,2))
# plot(1:k,fdp_seq,type="b",col="red",ylim=c(0,1),pch = 16,
#      xlab = "Threshold c",ylab = "Empirical FDR/Power")
# lines(1:k,pw_seq,type="b",col="blue",ylim=c(0,1),pch = 16)
# abline(h = q, lty = 2)
# legend("topright",legend=c("Empirical power","Empirical FDP"),pch=c(16,16),
#        col=c("blue","red"),lwd=2,lty=c(1,1))
# 
# size = result_adp$size[1:7]
# ratio = result_adp$ratio[1:7]
# 
# plot(size,col="blue",type="b",ylim=c(0,35),xlab="Threshold c",ylab="Subset size and 20*log(Ratio)",pch=15)
# lines(20*rt,type="b",col="red",pch=15)
# abline(h=20,lty=2)
# abline(v=4,lty=2,col="green")
# legend("topright",legend=c("Subset size","20*log(Ratio)"),pch=c(15,15),
#        col=c("blue","red"),lwd=2,lty=c(1,1))






#study different k: number of machines
repN = 10
d = 50
kset = c(1,2,5,8,10,20)
sp = 20

detectCores()
cl<- makeCluster(detectCores())
registerDoParallel(cl)
res = NULL
for(ik in kset){
    ires = foreach(rn = 1:repN, .combine="+") %dopar% simu_i(d,ik,sp)
    res = rbind(res, ires/repN)
}
stopCluster(cl)
FDP = res[,1]; PW = res[,(7+1)]; FDP_min = res[,2]; PW_min = res[,(7+2)];
ifdp = res[,3]; ipw = res[,(7+3)];
ufdp = res[,4]; upw = res[,(7+4)]; mfdp = res[,5]; mpw = res[,(7+5)];
Xfdp = res[,6]; Xpw = res[,(7+6)]; Afdp = res[,7]; Apw = res[,(7+7)]

plot(kset, FDP, col = "red",pch = 16,
     xlab = "k: number of machines",
     ylab = "Empirical FDR",
     ylim = c(0,1), type = "b")
lines(kset, ifdp, col = "blue",pch = 16, type = "b")
lines(kset, ufdp, col = "green", pch = 16, type = "b")
lines(kset, mfdp, col = "brown", pch = 16, type = "b")
lines(kset, Xfdp, col = "orange", pch = 15, type = "b")
lines(kset, Afdp, col = "pink", pch = 15, type = "b")
lines(kset, FDP_min, col = "red", pch = 15, type = "b", lty = 2)
abline(h = q, lty = 2)
legend("topleft",legend=c("Adaptive_1","Adaptive_2","Intersection","Union",
                           "Median","Xie","AKO"),
       pch=c(16,15,16,16,16,15,15),
       col=c("red","red", "blue","green",
             "brown","orange","pink"),lwd=1,lty=c(1,2,1,1,1,1,1))

plot(kset, PW, col = "red",pch = 16,
     xlab = "k: number of machines",
     ylab = "Empirical power",
     ylim = c(0,1), type = "b")
lines(kset, ipw, col = "blue",pch = 16, type = "b")
lines(kset, upw, col = "green", pch = 16, type = "b")
lines(kset, mpw, col = "brown", pch = 16, type = "b")
lines(kset, Xpw, col = "orange", pch = 15, type = "b")
lines(kset, Apw, col = "pink", pch = 15, type = "b")
lines(kset, PW_min, col = "red", pch = 15, type = "b", lty=2)

##study different d: dimensions
repN = 100
n = 1000
k = 10
sp = 10


dset = c(15,30,45,60,75,90)
ld = length(dset)
detectCores()
cl<- makeCluster(detectCores())
registerDoParallel(cl)
res = NULL
for(id in dset){
    ires = foreach(rn = 1:repN, .combine="+") %dopar% simu_i(id,k,sp)
    res = rbind(res, ires/repN)
}
stopCluster(cl)
FDP = res[,1]; PW = res[,(7+1)]; FDP_min = res[,2]; PW_min = res[,(7+2)];
ifdp = res[,3]; ipw = res[,(7+3)];
ufdp = res[,4]; upw = res[,(7+4)]; mfdp = res[,5]; mpw = res[,(7+5)];
Xfdp = res[,6]; Xpw = res[,(7+6)]; Afdp = res[,7]; Apw = res[,(7+7)]



par(mfrow=c(1,2))
plot(dset, FDP, col = "red",pch = 16,
     xlab = "d: dimension", xlim = c(10,90),
     ylab = "Empirical FDR",
     ylim = c(0,1), type = "b")
lines(dset, ifdp, col = "blue",pch = 16, type = "b")
lines(dset, ufdp, col = "green", pch = 16, type = "b")
lines(dset, mfdp, col = "brown", pch = 16, type = "b")
lines(dset, Xfdp, col = "orange", pch = 15, type = "b")
lines(dset, Afdp, col = "pink", pch = 15, type = "b")
lines(dset, FDP_min, col = "red", pch = 15, type = "b", lty = 2)
abline(h = q, lty = 2)
legend("topleft",legend=c("ADAGES","ADAGES_m","Intersection","Union",
                           "Median","Xie","AKO"),
       pch=c(16,15,16,16,16,15,15),
       col=c("red","red", "blue","green",
             "brown","orange","pink"),lwd=1,lty=c(1,2,1,1,1,1,1),
       ncol = 1)

plot(dset, PW, col = "red",pch = 16,
     xlab = "d: dimension", xlim = c(10,90),
     ylab = "Empirical power",
     ylim = c(0,1), type = "b")
lines(dset, ipw, col = "blue",pch = 16, type = "b")
lines(dset, upw, col = "green", pch = 16, type = "b")
lines(dset, mpw, col = "brown", pch = 16, type = "b")
lines(dset, Xpw, col = "orange", pch = 15, type = "b")
lines(dset, Apw, col = "pink", pch = 15, type = "b")
lines(dset, PW_min, col = "red", pch = 15, type = "b",lty=2)




#machine-wise
# repN = 1
# n = 1000
# k = 10
# sp = 10
# 
# 
# d = 80
# 
# count = matrix(rep(0,repN*(k+1)),nrow = repN)
# 
# id = d
# for(rp in 1:repN){
#     #distributed generation
#     #sp = 0.4*id
#     nonzero = sample(id,sp)
#     beta = numeric(id)
#     indv_set = NULL
#     indv_set1 = NULL
#     p_val = NULL
#     nn = n/k
#     beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
#     for(i in 1:k){
#         mu = rep(0,id)
#         rho = 0.25
#         Sigma = toeplitz(rho^(0:(id-1)))
#         X = matrix(rnorm(nn*id),nn) %*% chol(Sigma)
#         y.sample = function(X) X %*% beta + rnorm(nn)
#         y = y.sample(X)
#         result = knockoff.filter(X, y,
#                                  knockoffs = create.second_order,
#                                  statistic = stat.glmnet_lambdasmax,
#                                  fdr = q)
#         selected_set = result$selected
#         indv_set[[i]] = result$selected
#         # i_fdp = sum(beta[selected_set] == 0) / max(1, length(selected))
#         # i_pw = sum(beta[selected_set] != 0) / length(nonzero)
#         # ifdp[rn,i] = i_fdp
#         # ipw[rn,i] = i_pw
#         result1 = knockoff.filter(X, y,
#                                   knockoffs = create.second_order,
#                                   statistic = stat.glmnet_lambdasmax,
#                                   fdr = q/k)
#         selected_set1 = result1$selected
#         indv_set1[[i]] = selected_set1
#         p_val[[i]] = p_values(result$statistic)
#     }
# 
#     #union
#     kk = 1
#     result = agg(k,indv_set,kk)
#     selected = result$set
#     fdp = sum(beta[selected] == 0) / max(1, length(selected))
#     pw = sum(beta[selected] != 0) / length(nonzero)
#     ufdp = fdp
#     upw = pw
#     m = result$m
#     for(ic in 0:k){
#         count[rp,ic+1] = sum(m==ic)
#     }
# }
# 
# boxplot(count,border = "blue",
#      xlab = "statistics m_j",
#      ylim = c(0,max(count)+5),
#      xaxt = "n")
# axis(1,1:(k+1),0:k)
# 
# 
# #result_addaptive
# result_adp = agg_adp(k,indv_set)
# c = result_adp$c
# cset = result_adp$cset
# selected = result_adp$set
# fdp = sum(beta[selected] == 0) / max(1, length(selected))
# pw = sum(beta[selected] != 0) / length(nonzero)
# #print(paste0("FDP = ",fdp,", Power = ",pw))
# FDP = fdp
# PW = pw
# 
# #mean
# kk = floor((k+1)/2)
# result = agg(k,indv_set,kk)
# selected = result$set
# fdp = sum(beta[selected] == 0) / max(1, length(selected))
# pw = sum(beta[selected] != 0) / length(nonzero)
# mfdp = fdp
# mpw = pw
# 
# #intersection
# kk = k
# result = agg(k,indv_set,kk)
# selected = result$set
# fdp = sum(beta[selected] == 0) / max(1, length(selected))
# pw = sum(beta[selected] != 0) / length(nonzero)
# ifdp = fdp
# ipw = pw
# 
# #Xie
# kk = 1
# result = agg(k,indv_set1,kk)
# selected = result$set
# fdp = sum(beta[selected] == 0) / max(1, length(selected))
# pw = sum(beta[selected] != 0) / length(nonzero)
# Xfdp = fdp
# Xpw = pw
# 
# #AKO
# AKO_select = AKO(p_val, d, k)
# fdp = sum(beta[AKO_select] == 0) / max(1, length(AKO_select))
# pw = sum(beta[AKO_select] != 0) / length(nonzero)
# Afdp = fdp
# Apw = pw
# 
# 
# mw_fdr = numeric(k)
# mw_pw = numeric(k)
# for(i in 1:k){
#     selected = indv_set[[i]]
#     fdp = sum(beta[selected] == 0) / max(1, length(selected))
#     pw = sum(beta[selected] != 0) / length(nonzero)
#     mw_fdr[i] = fdp
#     mw_pw[i] = pw
# }
# 
# par(mfrow = c(1,2))
# barplot(c(mw_fdr,FDP,ifdp,ufdp,mfdp,Xfdp,Afdp),
#         col = c(rep("gray",k),
#                 "red","blue","green","brown","orange","pink"),
#         ylim = c(0,1))
# abline(h = q, lty = 2)
# legend("topleft",legend=c(paste0("m1-m",k),
#          "Adaptive","Intersection","Union","Median","Xie","AKO"),
#        col = c("gray",
#                "red","blue","green","brown","orange","pink"),
#        pch = rep(15,7),
#        ncol = 2)
# 
# barplot(c(mw_pw,PW,ipw,upw,mpw,Xpw,Apw),
#         col = c(rep("gray",k),
#                 "red","blue","green","brown","orange","pink"),
#         ylim = c(0,1))
# 
# 
# 

# freq
# comp_fq = function(d,k,sp,nonzero,j0){
#     library(knockoff)
#     beta = numeric(d) 
#     indv_set = NULL
#     p_val = NULL
#     indv_set1 = NULL
#     nn = n/k
#     beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
#     for(i in 1:k){
#         mu = rep(0,d)
#         rho = 0.25
#         Sigma = toeplitz(rho^(0:(d-1)))
#         X = matrix(rnorm(nn*d),nn) %*% chol(Sigma)
#         y.sample = function(X) X %*% beta + rnorm(nn)
#         y = y.sample(X)
#         result = knockoff.filter(X, y, 
#                                  knockoffs = create.second_order, 
#                                  statistic = stat.glmnet_lambdasmax,
#                                  fdr = q)
#         selected_set = result$selected
#         indv_set[[i]] = result$selected
#         p_val[[i]] = p_values(result$statistic)
#     }
#     result_adp = agg_adp(k,indv_set)
#     c = result_adp$c
#     result_min = agg_min(k,indv_set)
#     c_min = result_min$c
#     f_seq = numeric(k)
#     for(i in 1:k) f_seq[i] = sum(indv_set[[i]]==j0)
#     return(c(sum(f_seq),c,c_min))
# }
# 
# 
# nonzero = sample(d,sp)
# j0 = nonzero[1]
# detectCores()
# cl<- makeCluster(detectCores())
# registerDoParallel(cl)
# d = 50
# k = 10
# sp = 10
# repN = 500
# frq = numeric(repN)
# res = foreach(rr = 1:repN, .combine = "rbind") %dopar% comp_fq(d,k,sp,nonzero,j0)
# stopCluster(cl)
# fq = res[,1]
# c = res[,2]
# c_min = res[,3]
# par(mfrow = c(1,2))
# hist(fq, 40, xlim = c(0,k), 
#      main = "Histogram of frequency: fixed j_0",
#      col = "orange", border = "orange",
#      ylab = "frequence", xlab = "selected by #machines: fixed nonzero j_0")
# boxplot(cbind(c,c_min), ylim = c(0,k), main = "Boxplot of adaptive thresholds",
#         border = c("blue","red"), ylab = "adaptive threshold")
# abline(h = 1, col = "orange", lty = 2)
# abline(h = k/2, col = "orange", lty = 2)
# 
