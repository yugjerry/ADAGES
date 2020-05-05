library(knockoff)

repN = 100

#model_parameters
n = 1000
d = 50
k = 10
sp = 20
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

FDP = numeric(repN)
PW = numeric(repN)
ifdp = matrix(rep(0,repN*k),nrow = repN)
ipw = matrix(rep(0,repN*k),nrow = repN)
fdp_seq = numeric(k)
pw_seq = numeric(k)
ssize = numeric(k)
for(rn in 1:repN){
    #distributed generation
    indv_set = NULL
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
                                 knockoffs = create.fixed, 
                                 statistic = stat.glmnet_lambdasmax,
                                 fdr = q)
        selected_set = result$selected
        indv_set[[i]] = result$selected
        i_fdp = sum(beta[selected_set] == 0) / max(1, length(selected))
        i_pw = sum(beta[selected_set] != 0) / length(nonzero)
        ifdp[rn,i] = i_fdp
        ipw[rn,i] = i_pw
    }
    
    #result_addaptive
    result_adp = agg_adp(k,indv_set)
    c = result_adp$c
    cset = result_adp$cset
    selected = result_adp$set
    fdp = sum(beta[selected] == 0) / max(1, length(selected))
    pw = sum(beta[selected] != 0) / length(nonzero)
    #print(paste0("FDP = ",fdp,", Power = ",pw))
    FDP[rn] = fdp
    PW[rn] = pw
    
    #fixed_threshold
    for(kk in 1:k){
        result = agg(k,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        fdp_seq[kk] = fdp + fdp_seq[kk]
        pw_seq[kk] = pw + pw_seq[kk]
        ssize[kk] = length(selected)
    }
    
}

fdp_seq = fdp_seq/repN
pw_seq = pw_seq/repN

fdpi = fdp_seq[k]
pwi = pw_seq[k]
fdpu = fdp_seq[1]
pwu = pw_seq[1]

fdp_adp = mean(FDP)
pw_adp = mean(PW)

par(mfrow=c(1,2))
plot(1:k,fdp_seq,type="b",col="red",ylim=c(0,1),pch = 16,
     xlab = "Threshold c",ylab = "Empirical FDR/Power")
lines(1:k,pw_seq,type="b",col="blue",ylim=c(0,1),pch = 16)
abline(h = q, lty = 2)
legend("topright",legend=c("Empirical power","Empirical FDP"),pch=c(16,16),
       col=c("blue","red"),lwd=2,lty=c(1,1))

size = result_adp$size[1:7]
ratio = result_adp$ratio[1:7]

plot(size,col="blue",type="b",ylim=c(0,35),xlab="Threshold c",ylab="Subset size and 20*log(Ratio)",pch=15)
lines(20*rt,type="b",col="red",pch=15)
abline(h=20,lty=2)
abline(v=4,lty=2,col="green")
legend("topright",legend=c("Subset size","20*log(Ratio)"),pch=c(15,15),
       col=c("blue","red"),lwd=2,lty=c(1,1))






##study different k: number of machines
repN = 100

kset = c(1,2,5,8,10,20)
lk = length(kset)

FDP = numeric(lk)
PW = numeric(lk)
ifdp = numeric(lk)
ipw = numeric(lk)
ufdp = numeric(lk)
upw = numeric(lk)
Xfdp = numeric(lk)
Xpw = numeric(lk)
mfdp = numeric(lk)
mpw = numeric(lk)

for(rn in 1:repN){
    set.seed(rn)
    for(iik in 1:lk){
        ik = kset[iik]
        #distributed generation
        
        indv_set = NULL
        indv_set1 = NULL
        nn = n/ik
        beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
        for(i in 1:ik){
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
                                     fdr = q/ik)
            selected_set1 = result1$selected
            indv_set1[[i]] = selected_set1
        }
        
        #result_addaptive
        result_adp = agg_adp(ik,indv_set)
        c = result_adp$c
        cset = result_adp$cset
        selected = result_adp$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        #print(paste0("FDP = ",fdp,", Power = ",pw))
        FDP[iik] = fdp + FDP[iik]
        PW[iik] = pw + PW[iik]
        
        #union
        kk = 1
        result = agg(ik,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        ufdp[iik] = fdp + ufdp[iik]
        upw[iik] = pw + upw[iik]
        
        #intersection
        kk = ik
        result = agg(ik,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        ifdp[iik] = fdp + ifdp[iik]
        ipw[iik] = pw + ipw[iik]
        
        #Median
        kk = floor((ik+1)/2)
        result = agg(ik,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        mfdp[iik] = fdp + mfdp[iik]
        mpw[iik] = pw + mpw[iik]
        
        #Xie
        kk = 1
        result = agg(ik,indv_set1,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        Xfdp[iik] = fdp + Xfdp[iik]
        Xpw[iik] = pw + Xpw[iik]
    }
}

FDP = FDP/repN
PW = PW/repN
ifdp = ifdp/repN
ipw = ipw/repN
ufdp = ufdp/repN
upw = upw/repN
Xfdp = Xfdp/repN
Xpw = Xpw/repN
mfdp = mfdp/repN
mpw = mpw/repN

plot(kset, FDP, col = "red",pch = 16, 
     xlab = "k: number of machines", 
     ylab = "Empirical FDR",
     ylim = c(0,1), type = "b")
lines(kset, ifdp, col = "blue",pch = 16, type = "b")
lines(kset, ufdp, col = "green", pch = 16, type = "b")
lines(kset, mfdp, col = "brown", pch = 16, type = "b")
lines(kset, Xfdp, col = "orange", pch = 15, type = "b")
abline(h = q, lty = 2)
legend("topright",legend=c("Adaptive","Intersection","Union","Median","Xie"),
       pch=c(16,16,16,16,15),
       col=c("red","blue","green","brown","orange"),lwd=1,lty=c(1,1,1,1,1))

plot(kset, PW, col = "red",pch = 16, 
     xlab = "k: number of machines", 
     ylab = "Empirical power",
     ylim = c(0,1), type = "b")
lines(kset, ipw, col = "blue",pch = 16, type = "b")
lines(kset, upw, col = "green", pch = 16, type = "b")
lines(kset, mpw, col = "brown", pch = 16, type = "b")
lines(kset, Xpw, col = "orange", pch = 15, type = "b")



##study different d: dimensions
repN = 1
n = 1000
k = 5
sp = 10


dset = c(20)
ld = length(dset)

FDP = numeric(ld)
PW = numeric(ld)
ifdp = numeric(ld)
ipw = numeric(ld)
ufdp = numeric(ld)
upw = numeric(ld)
Xfdp = numeric(ld)
Xpw = numeric(ld)
mfdp = numeric(ld)
mpw = numeric(ld)

for(rn in 1:repN){
    for(iid in 1:ld){
        id = dset[iid]
        #distributed generation
        #sp = 0.4*id
        nonzero = sample(id,sp)
        beta = numeric(id)
        indv_set = NULL
        indv_set1 = NULL
        nn = n/k
        beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
        for(i in 1:k){
            mu = rep(0,id)
            rho = 0.25
            Sigma = toeplitz(rho^(0:(id-1)))
            X = matrix(rnorm(nn*id),nn) %*% chol(Sigma)
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
        }
        
        #result_addaptive
        result_adp = agg_adp(k,indv_set)
        c = result_adp$c
        cset = result_adp$cset
        selected = result_adp$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        #print(paste0("FDP = ",fdp,", Power = ",pw))
        FDP[iid] = fdp + FDP[iid]
        PW[iid] = pw + PW[iid]
        
        #union
        kk = 1
        result = agg(k,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        ufdp[iid] = fdp + ufdp[iid]
        upw[iid] = pw + upw[iid]
        
        #mean
        kk = floor((k+1)/2)
        result = agg(k,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        mfdp[iid] = fdp + mfdp[iid]
        mpw[iid] = pw + mpw[iid]
        
        #intersection
        kk = k
        result = agg(k,indv_set,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        ifdp[iid] = fdp + ifdp[iid]
        ipw[iid] = pw + ipw[iid]
        
        #Xie
        kk = 1
        result = agg(k,indv_set1,kk)
        selected = result$set
        fdp = sum(beta[selected] == 0) / max(1, length(selected))
        pw = sum(beta[selected] != 0) / length(nonzero)
        Xfdp[iid] = fdp + Xfdp[iid]
        Xpw[iid] = pw + Xpw[iid]
    }
}

FDP = FDP/repN
PW = PW/repN
ifdp = ifdp/repN
ipw = ipw/repN
ufdp = ufdp/repN
upw = upw/repN
Xfdp = Xfdp/repN
Xpw = Xpw/repN
mfdp = mfdp/repN
mpw = mpw/repN

par(mfrow=c(1,2))
plot(dset, FDP, col = "red",pch = 16, 
     xlab = "d: dimension", xlim = c(10,90),
     ylab = "Empirical FDR",
     ylim = c(0,1), type = "b")
lines(dset, ifdp, col = "blue",pch = 16, type = "b")
lines(dset, ufdp, col = "green", pch = 16, type = "b")
lines(dset, mfdp, col = "brown", pch = 16, type = "b")
lines(dset, Xfdp, col = "orange", pch = 15, type = "b")
abline(h = q, lty = 2)
legend("topright",legend=c("Adaptive","Intersection","Union","Median","Xie"),
       pch=c(16,16,16,16,15),
       col=c("red","blue","green","brown","orange"),lwd=1,lty=c(1,1,1,1,1))

plot(dset, PW, col = "red",pch = 16, 
     xlab = "d: dimension", xlim = c(10,90),
     ylab = "Empirical power",
     ylim = c(0,1), type = "b")
lines(dset, ipw, col = "blue",pch = 16, type = "b")
lines(dset, upw, col = "green", pch = 16, type = "b")
lines(dset, mpw, col = "brown", pch = 16, type = "b")
lines(dset, Xpw, col = "orange", pch = 15, type = "b")




#machine-wise
repN = 100
n = 1000
k = 10
sp = 10


d = 60

count = matrix(rep(0,repN*(k+1)),nrow = repN)

id = d
for(rp in 1:repN){
    #distributed generation
    #sp = 0.4*id
    nonzero = sample(id,sp)
    beta = numeric(id)
    indv_set = NULL
    indv_set1 = NULL
    nn = n/k
    beta[nonzero] = 2 * (2*rbinom(sp,1,0.5) - 1)
    for(i in 1:k){
        mu = rep(0,id)
        rho = 0.25
        Sigma = toeplitz(rho^(0:(id-1)))
        X = matrix(rnorm(nn*id),nn) %*% chol(Sigma)
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
    }
    
    #union
    kk = 1
    result = agg(k,indv_set,kk)
    selected = result$set
    fdp = sum(beta[selected] == 0) / max(1, length(selected))
    pw = sum(beta[selected] != 0) / length(nonzero)
    ufdp = fdp 
    upw = pw
    m = result$m
    for(ic in 0:k){
        count[rp,ic+1] = sum(m==ic)
    }
}

boxplot(count,border = "blue",
     xlab = "statistics m_j",
     ylim = c(0,max(count)+5),
     xaxt = "n")
axis(1,1:(k+1),0:k)


#result_addaptive
result_adp = agg_adp(k,indv_set)
c = result_adp$c
cset = result_adp$cset
selected = result_adp$set
fdp = sum(beta[selected] == 0) / max(1, length(selected))
pw = sum(beta[selected] != 0) / length(nonzero)
#print(paste0("FDP = ",fdp,", Power = ",pw))
FDP = fdp
PW = pw 

#mean
kk = floor((k+1)/2)
result = agg(k,indv_set,kk)
selected = result$set
fdp = sum(beta[selected] == 0) / max(1, length(selected))
pw = sum(beta[selected] != 0) / length(nonzero)
mfdp = fdp 
mpw = pw 

#intersection
kk = k
result = agg(k,indv_set,kk)
selected = result$set
fdp = sum(beta[selected] == 0) / max(1, length(selected))
pw = sum(beta[selected] != 0) / length(nonzero)
ifdp = fdp 
ipw = pw

#Xie
# kk = 1
# result = agg(k,indv_set1,kk)
# selected = result$set
# fdp = sum(beta[selected] == 0) / max(1, length(selected))
# pw = sum(beta[selected] != 0) / length(nonzero)
# Xfdp = fdp
# Xpw = pw 


mw_fdr = numeric(k)
mw_pw = numeric(k)
for(i in 1:k){
    selected = indv_set[[i]]
    fdp = sum(beta[selected] == 0) / max(1, length(selected))
    pw = sum(beta[selected] != 0) / length(nonzero)
    mw_fdr[i] = fdp
    mw_pw[i] = pw
}

par(mfrow = c(1,2))
barplot(c(mw_fdr,FDP,ifdp,ufdp,mfdp,Xfdp),
        col = c(rep("gray",k),
                "red","blue","green","brown","orange"),
        ylim = c(0,1))
abline(h = q, lty = 2)
legend("topleft",legend=matrix(c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10",
         "Adaptive","Intersection","Union","Median","Xie"),nrow=5),
       col = matrix(c(rep("gray",k),
               "red","blue","green","brown","orange"),nrow=5),
       pch = matrix(rep(15,10),nrow=5),
       ncol = 3)

barplot(c(mw_pw,PW,ipw,upw,mpw,Xpw),
        col = c(rep("gray",k),
                "red","blue","green","brown","orange"),
        ylim = c(0,1))



