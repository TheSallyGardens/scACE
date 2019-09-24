## EM implementation for "Jointly clustering chromatin accessibility and gene expression data with application to single cell genomics"

## calculate prob_z_atac
calculate_prob_z_atac<-function(phi_atac,f1,f0,w_acc,qi){
  temp<-array(rep(0,i0*k0*p0),dim=c(i0,k0,p0))
  for(k in c(1:k0)){
    temp[,k,]<-t(t(f1*qi+f0*(1-qi))*(w_acc[k,]))+t(t(f0)*(1-w_acc[k,]))
  }
  prob_z_atac<-t(t(apply(log(temp),c(1,2),sum))+log(phi_atac))
  prob_z_atac2<-prob_z_atac
  for(k in c(1:k0)){
    if (k0<=2){
      prob_z_atac[,k] <- 1/(1+(exp(prob_z_atac2[,-k]-prob_z_atac2[,k])))
    }
    else{
      prob_z_atac[,k] <- 1/(1+(rowSums(exp(prob_z_atac2[,-k]-prob_z_atac2[,k]))))
    }
  }
  return(prob_z_atac)
}

## calculate prob_z_u_atac
calculate_prob_z_u_atac<-function(prob_z_atac,f1,f0,w_acc,qi){
  temp<-array(c(1:(k0*i0*p0)),dim=c(i0,k0,p0))
  temp2<-array(c(1:(k0*i0*p0)),dim=c(i0,k0,p0))
  for(k in c(1:k0)){
    temp[,k,]<-t((t(f1*qi)+t(f0*(1-qi)))*(w_acc[k,]))
    temp2[,k,]<-t(t(f0)*(1-w_acc[k,]))
  }
  prob_z_u_atac<-array(rep(prob_z_atac,p0),dim=c(i0,k0,p0))/(1+exp(log(temp2)-log(temp)))
  return(prob_z_u_atac)
}

## calculate prob_z_u_ut_atac
calculate_prob_z_u_ut_atac<-function(prob_z_u_atac,f1,f0,qi){
  prob_z_u_ut_atac<-prob_z_u_atac
  for(k in c(1:k0)){
    prob_z_u_ut_atac[,k,]<-prob_z_u_atac[,k,]/(1+exp(log(f0*(1-qi))-log(f1*qi)))
  }
  return(prob_z_u_ut_atac)
}

## calculate prob_z_rna
calculate_prob_z_rna<-function(phi_rna,g1,g0,w_exp,pi_exp,ql,po){
  temp<-array(rep(0,l0*k0*p0),dim=c(l0,k0,p0))
  for(k in c(1:k0)){
    temp[,k,1:po]<-t(t(g1[,1:po]*pi_exp[1,]*ql)*w_exp[k,1:po]+t(g1[,1:po]*pi_exp[2,]*ql)*(1-w_exp[k,1:po])+t(g0[,1:po]*(1-pi_exp[1,]))*w_exp[k,1:po]+t(g0[,1:po]*(1-pi_exp[2,]))*(1-w_exp[k,1:po])+t(pi_exp[1,]*g0[,1:po]*(1-ql))*w_exp[k,1:po]+(1-w_exp[k,1:po])*t(pi_exp[2,]*g0[,1:po]*(1-ql)))
    temp[,k,(po+1):p0]<-t(t(ql*g1[,(po+1):p0])*w_exp[k,(po+1):p0]+t((1-ql)*g0[,(po+1):p0])*w_exp[k,(po+1):p0]+t(g0[,(po+1):p0])*(1-w_exp[k,(po+1):p0]))
  }
  prob_z_rna<-t(t(apply(log(temp),c(1,2),sum))+log(phi_rna))
  prob_z_rna2<-prob_z_rna
  for(k in c(1:k0)){
    if (k0<=2){
      prob_z_rna[,k]<-1/(1+(exp(prob_z_rna2[,-k]-prob_z_rna2[,k])))
    }
    else{
      prob_z_rna[,k]<-1/(1+(rowSums(exp(prob_z_rna2[,-k]-prob_z_rna2[,k]))))
    }
  }

  return(prob_z_rna)
}

## calculate prob_z_u_rna
calculate_prob_z_u_rna<-function(prob_z_rna,g1,g0,w_exp,pi_exp,ql,po){
  prob_z_u_rna<-array(rep(0,(l0*k0*po)),dim=c(l0,k0,po))
  temp<-array(rep(0,(l0*k0*po)),dim=c(l0,k0,po))
  temp2<-array(rep(0,(l0*k0*po)),dim=c(l0,k0,po))
  for(k in c(1:k0)){
    temp[,k,]<-t(log(w_exp[k,1:po])+t(log(pi_exp[1,]*ql*g1[,1:po]+pi_exp[1,]*(1-ql)*g0[,1:po]+g0[,1:po]*(1-pi_exp[1,]))))
    temp2[,k,]<-t(log(1-w_exp[k,1:po])+t(log(pi_exp[2,]*ql*g1[,1:po]+pi_exp[2,]*(1-ql)*g0[,1:po]+g0[,1:po]*(1-pi_exp[2,]))))
    prob_z_u_rna[,k,]<-matrix(rep(prob_z_rna[,k],po),l0,po)/(1+exp(temp2[,k,]-temp[,k,]))
  }
  return(prob_z_u_rna)
}

## calculate prob_z_u_v_rna
calculate_prob_z_u_v_rna<-function(prob_z_u_rna,g1,g0,pi_exp,ql,po){
  prob_z_u_v_rna<-array(rep(0,(l0*k0*po)),dim=c(l0,k0,po))
  deno<-(1+exp(log((1-pi_exp[1,])*g0[,1:po])-log(pi_exp[1,]*(g0[,1:po]*(1-ql)+g1[,1:po]*ql))))
  for(k in c(1:k0)){
    prob_z_u_v_rna[,k,]<-prob_z_u_rna[,k,]/deno
  }
  return(prob_z_u_v_rna)
}

## calculate prob_z_1_u_v
calculate_prob_z_1_u_v<-function(prob_z_rna,prob_z_u_rna,g1,g0,pi_exp,ql,po){
  prob_z_1_u_v<-array(rep(0,(l0*k0*po)),dim=c(l0,k0,po))
  deno<-(1+exp(log((1-pi_exp[2,])*g0[,1:po])-log(pi_exp[2,]*(g0[,1:po]*(1-ql)+g1[,1:po]*ql))))
  for(k in c(1:k0)){
    prob_z_1_u_v[,k,]<-(prob_z_rna[,k]-prob_z_u_rna[,k,])/deno
  }
  return(prob_z_1_u_v)
}

## calculate prob_z_v
calculate_prob_z_v<-function(prob_z_rna,prob_z_u_v,prob_z_1_u_v,w_exp,g1,g0,ql,po){
  prob_z_v<-array(rep(0,l0*k0*p0),dim=c(l0,k0,p0))
  prob_z_v[,,1:po]<-prob_z_u_v+prob_z_1_u_v
  for(k in c(1:k0)){
    prob_z_v[,k,(po+1):p0]<-matrix(rep(prob_z_rna[,k],p0-po),l0,p0-po)/(1+exp(t(log((1-w_exp[k,(po+1):p0])*t(g0[,(po+1):p0])))-t(log((w_exp[k,(po+1):p0]*t(ql*g1[,(po+1):p0]+(1-ql)*g0[,(po+1):p0]))))))
  }
  return(prob_z_v)
}

## calculate prob_z_v_vt
calculate_prob_z_v_vt_rna<-function(prob_z_v,g1,g0,ql){
  prob_z_v_vt<-prob_z_v
  mult<-(ql*g1)/(ql*g1+(1-ql)*g0)
  for(k in c(1:k0)){
    prob_z_v_vt[,k,]<-prob_z_v[,k,]*mult
  }
  return(prob_z_v_vt)
}

## sum function across different dimension
sum_k_g<-function(prob,po=p0){
  return(apply(prob[,,1:po],1,sum))
}

sum_k_gi<-function(prob,po=p0){
  return(apply(prob[,,1:po],1,sum))

}

sum_prob_z<-function(prob_z){
  return(apply(prob_z,2,sum))
}

sum_prob_zi<-function(prob_z){
  return(apply(prob_z,2,sum))
}

sum_prob_z_u<-function(prob_z_u,po){
  return(apply(prob_z_u[,,1:po],c(2,3),sum))
}

sum_prob_z_ui<-function(prob_z_u){
  return(apply(prob_z_u,c(2,3),sum))
}

## update pi_exp
update_pi3<-function(sum_prob_z_u_v_rna,sum_prob_z_1_u_v_rna,sum_prob_z_u_rna,pi_exp,po){
  temp<-1/(1+exp(log(sum_prob_z_u_rna-sum_prob_z_u_v_rna)-log(sum_prob_z_u_v_rna)))
  temp[which(temp>=1-10^-5)]=1-10^-5

  tmp<-1/(1-1/sum_prob_z_1_u_v_rna+exp(log(po-sum_prob_z_u_rna-sum_prob_z_1_u_v_rna)-log(sum_prob_z_1_u_v_rna)))

  flag<-which(tmp<=temp)
  pi_exp[1,flag]<-temp[flag]
  pi_exp[2,flag]<-tmp[flag]
  return(pi_exp)
}

## update phi_atac
update_phi_atac<-function(phi_atac,f1,f0,w_acc,qi){
  prob_z_atac<-calculate_prob_z_atac(phi_atac,f1,f0,w_acc,qi)
  temp<-sum_prob_zi(prob_z_atac)
  phi_atac<-(1+temp)/(i0+k0)
  return(phi_atac)
}

## update phi_rna
update_phi_rna<-function(phi_rna,g1,g0,w_exp,pi_exp,ql,po){
  prob_z_rna<-calculate_prob_z_rna(phi_rna,g1,g0,w_exp,pi_exp,ql,po)
  temp<-sum_prob_z(prob_z_rna)
  phi_rna<-(1+temp)/(l0+k0)
  return(phi_rna)
}

## update w_exp
update_w_exp<-function(sum_prob_z_rna,sum_prob_z_u_rna,sum_prob_z_v_rna,w_acc,w_exp,phi_1,po){
  w_exp[,1:po]<-(1+sign(w_acc[,1:po]*phi_1-1)*exp(log(abs(w_acc[,1:po]*phi_1-1))-log(sum_prob_z_u_rna[,1:po])))/(exp(log(sum_prob_z_rna)-log(sum_prob_z_u_rna[,1:po]))+sign(phi_1-1)*exp(log(abs(phi_1-2))-log(sum_prob_z_u_rna[,1:po])))
  w_exp[which(w_exp>=1-10^-6)]<-1-10^-6
  w_exp[which(w_exp<=1-10^-6)]<-10^-6
  w_exp[,(po+1):p0]<-(1+sum_prob_z_v_rna[,(po+1):p0])/(2+sum_prob_z_rna)
  return(w_exp)
}

## update w_acc
update_w_acc_grid_search_outer<-function(sum_prob_z_acc,sum_prob_z_u_acc,w_acc,w_exp,phi_1,po){
  # w_acc seq
  temp <- seq(from=0.01, to=1, by=0.01)
  ntemp <- length(temp)

  w_acc_temp <- matrix(rep(temp,each=k0*po),nrow=k0*po,byrow = FALSE)
  w_exp_linked <- as.vector(w_exp[,(1:po)])
  sum_prob_z_u_acc_linked <- as.vector(sum_prob_z_u_acc[,(1:po)])
  diff_z_and_z_u_acc <- as.vector(sum_prob_z_acc - sum_prob_z_u_acc[,(1:po)])

  pt1 <- outer(sum_prob_z_u_acc_linked,log(temp))
  pt2 <- outer(diff_z_and_z_u_acc,log(1-temp))
  pt3 <- outer(log(w_exp_linked),(temp*phi_1-1))
  pt4 <- outer(log(1-w_exp_linked),(phi_1-1-phi_1*temp))
  pt5 <- beta(w_acc_temp*phi_1,(phi_1-phi_1*w_acc_temp))

  y <- -pt1-pt2-pt3-pt4+log(pt5)
  w_acc_index <- apply(y,1, function(x) which.min(x))
  w_acc_all <- temp[w_acc_index]
  w_acc[,(1:po)] <- matrix(w_acc_all, nrow=k0, byrow = FALSE)

  for(k in c(1:k0)){
    w_acc[k,(po+1):p0]<-(1+1/sum_prob_z_u_acc[k,(po+1):p0])/(1+2/sum_prob_z_u_acc[k,(po+1):p0]+exp(log(sum_prob_z_acc[k]-sum_prob_z_u_acc[k,(po+1):p0])-log(sum_prob_z_u_acc[k,(po+1):p0])))
  }
  return(w_acc)
}

## update ql
update_ql<-function(sum_prob_z_v_vt,sum_prob_z_v){
  ql<-(1+1/sum_prob_z_v_vt)/(1+2/sum_prob_z_v_vt+exp(log(sum_prob_z_v-sum_prob_z_v_vt)-log(sum_prob_z_v_vt)))

  return(ql)
}

## update qi
update_qi<-function(sum_prob_z_u_ut,sum_prob_z_u_atac){

  qi<-(1+1/sum_prob_z_u_ut)/(1+2/sum_prob_z_u_ut+exp(log(sum_prob_z_u_atac-sum_prob_z_u_ut)-log(sum_prob_z_u_ut)))

  return(qi)
}

## update phi_1
update_phi_1_outer<-function(w_acc,w_exp,po){
  temp<-seq(0.2,50,0.2)
  w_acc0 <- as.vector(w_acc[,(1:po)])
  w_exp0 <- as.vector(w_exp[,(1:po)])
  prod1 <- outer(w_acc0,temp)
  prod2 <- outer((1-w_acc0),temp)

  pt1 <- (prod1-1)*log(w_exp0)
  pt2 <- (prod2-1)*log(1-w_exp0)
  pt3 <- beta(prod1,prod2)
  y <- pt1+pt2-log(pt3)
  res <- colSums(y)
  phi_1<-temp[which.max(res)]
  return(phi_1)
}

## update all the parameters
update_all2<-function(phi_atac,f1,f0,w_acc,phi_rna,g1,g0,w_exp,pi_exp,qi,ql,phi_1,po){
  prob_z_rna<-calculate_prob_z_rna(phi_rna,g1,g0,w_exp,pi_exp,ql,po)
  prob_z_u_rna<-calculate_prob_z_u_rna(prob_z_rna,g1,g0,w_exp,pi_exp,ql,po)
  prob_z_u_v<-calculate_prob_z_u_v_rna(prob_z_u_rna,g1,g0,pi_exp,ql,po)
  prob_z_1_u_v<-calculate_prob_z_1_u_v(prob_z_rna,prob_z_u_rna,g1,g0,pi_exp,ql,po)
  prob_z_v<-calculate_prob_z_v(prob_z_rna,prob_z_u_v,prob_z_1_u_v,w_exp,g1,g0,ql,po)
  prob_z_v_vt<-calculate_prob_z_v_vt_rna(prob_z_v,g1,g0,ql)

  phi_rna_new<-update_phi_rna(phi_rna,g1,g0,w_exp,pi_exp,ql,po)
  ql<-update_ql(sum_k_g(prob_z_v_vt),sum_k_g(prob_z_v))
  w_exp_new<-update_w_exp(sum_prob_z(prob_z_rna),sum_prob_z_u(prob_z_u_rna,po),sum_prob_z_u(prob_z_v,p0),w_acc,w_exp,phi_1,po)
  pi_exp_new<-update_pi3(sum_k_g(prob_z_u_v,po),sum_k_g(prob_z_1_u_v,po),sum_k_g(prob_z_u_rna,po),pi_exp,po)

  prob_z_atac<-calculate_prob_z_atac(phi_atac,f1,f0,w_acc,qi)
  prob_z_u_atac<-calculate_prob_z_u_atac(prob_z_atac,f1,f0,w_acc,qi)
  prob_z_u_ut<-calculate_prob_z_u_ut_atac(prob_z_u_atac,f1,f0,qi)

  phi_atac_new<-update_phi_atac(phi_atac,f1,f0,w_acc,qi)
  qi <- update_qi(sum_k_gi(prob_z_u_ut),sum_k_gi(prob_z_u_atac))
  w_acc <- update_w_acc_grid_search_outer(sum_prob_zi(prob_z_atac),sum_prob_z_ui(prob_z_u_atac),w_acc,w_exp_new,phi_1,po)

  phi_1<-update_phi_1_outer(w_acc,w_exp_new,po)

  # calculate posterior probability pst
  t1<-array(rep(0,i0*k0*p0),dim=c(i0,k0,p0))
  for(k in c(1:k0)){
    t1[,k,]<-t(t(qi*f1+(1-qi)*f0)*w_acc[k,]+t(f0)*(1-w_acc[k,]))
  }
  product<-apply(log(t1),c(1,2),sum)
  part1<-sum(product[,1])+sum(log(rowSums(t(t(exp(product-product[,1]))*phi_atac_new))))

  t2<-array(rep(0,l0*k0*p0),dim=c(l0,k0,p0))
  for(k in c(1:k0)){
    t2[,k,1:po]<-t( t(((g1[,1:po]*ql+g0[,1:po]*(1-ql))*pi_exp_new[1,])+(1-pi_exp_new[1,])*g0[,1:po])*w_exp_new[k,1:po]+t(((g1[,1:po]*ql+g0[,1:po]*(1-ql))*pi_exp_new[2,])+(1-pi_exp_new[2,])*g0[,1:po])*(1-w_exp_new[k,1:po]) )
    t2[,k,(po+1):p0]<-t( w_exp_new[k,(po+1):p0]*t(ql*g1[,(po+1):p0]+(1-ql)*g0[,(po+1):p0])+(1-w_exp_new[k,(po+1):p0])*t(g0[,(po+1):p0]) )
  }
  product2<-apply(log(t2),c(1,2),sum)
  a<-max(product2-product2[,1])
  part2<-sum(product2[,1])+a/2*l0+sum(log(rowSums(t(phi_rna_new*t(exp(product2-product2[,1]-a/2))))))


  part3<-sum((phi_1*w_acc[,1:po]-1)*log(w_exp_new[,1:po])+(phi_1-1-phi_1*w_acc[,1:po])*log(1-w_exp_new[,1:po])-log(beta(phi_1*w_acc[,1:po],phi_1-phi_1*w_acc[,1:po])))
  part3<-part3-sum(log(1-pi_exp_new[2,]))+sum(log(ql))+sum(log(1-ql))+sum(log(qi))+sum(log((1-qi)))

  for(k in c(1:k0)){
    part3<-part3+log(phi_atac_new[k])+log(phi_rna_new[k])
    part3<-part3+sum(log(w_acc[k,(po+1):p0])+log(1-w_acc[k,(po+1):p0])+log(w_exp_new[k,(po+1):p0])+log(1-w_exp_new[k,(po+1):p0]))
  }

  pst<-part1+part2+part3
  return(list(phi_atac=phi_atac_new,phi_rna=phi_rna_new,w_exp=w_exp_new,w_acc=w_acc,pi_exp=pi_exp_new,qi=qi,ql=ql,phi_1=phi_1,postprob=pst))
}

