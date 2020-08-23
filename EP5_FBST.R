#George Othon NUSP 103xxxxx
#Felipe Zaffalon NUSP 103xxxxx


n=20

#Tabela Hardy Weinberg 
x1=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,5,5,5,5,5,5,5,5,5,9,9,9,9,9,9,9,9)
x3=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7)



# Otimizando a função f

f = function (theta_1,i) {
  x2=20-x1[i]-x3[i]
  theta_3=(1-sqrt(theta_1))^2
  theta_2=1-theta_1-theta_3
  y=(theta_1^x1[i])*(theta_2^x2)*(theta_3^x3[i])
  return(y)
}




fs=rep(0, 36)
thetas=rep(0, 36)
i=1

while (i<=36){
  fmax=optimise(f, c(0,1),  maximum = TRUE, i)$object
  fs[i]=fmax
  theta_max=optimise(f, c(0,1), maximum = TRUE,i)$maximum
  thetas[i]=theta_max
  i=i+1
} 



# Hipotese não nula
g  = function(theta_1, theta_3, i){
  theta_2 = 1-theta_1-theta_3
  x2 = (20-x1[i]-x3[i])
  y=theta_1^x1[i]*theta_2^x2*theta_3^x3[i]
  return(y)
}

#MCMC com metropolis

mcmc = function(i){
  vetor=rep(0, 10000)
  n=2
  theta_1=1/3
  theta_3=1/3
  vetor[1]=c(g(theta_1,theta_3,i))
  
  while(1){
    pro_theta_1 = theta_1+rnorm(1,0,1)
    pro_theta_3 = theta_3+rnorm(1,0,1)
    pro_theta_2 = 1-pro_theta_1 - pro_theta_3
    
    if (pro_theta_1+pro_theta_3+pro_theta_2==1 & pro_theta_2>0  & pro_theta_1>0 & pro_theta_3>0){ #veifica se o ponto atende as condições
      alpha = min(1, g(pro_theta_1, pro_theta_3, i)/g(theta_1,theta_3, i)) #passo de Metropolis
      
      if (runif(1)< alpha){# Condição de aceitação
        vetor[n]= g(pro_theta_1,pro_theta_3, i)
        theta_1=pro_theta_1
        theta_3=pro_theta_3
      }
      else { # caso não for aceito
        vetor[n]= g(theta_1,theta_3, i)
      }
      n=n+1
    }
    if (n==10000){break}
  }
  return(vetor)
}

tabela_resultados=c()
i=1
while (i<=36){ #Rodar para todos os x1, x3
  t=1
  media=c()
  while (t<=5) { #Rodar 5 vezes, afim de buscar mais consistencia 
    f_max=fs[i]
    final = rep(0, 10000)
    atual = mcmc(i)
    j=1
    while (j<=length(agr)) {
      if (atual[j]>f_max){
        final[j]=1
      }
      else {
        final[j]=0
      }
      j=j+1
    }
    media[t]=1-mean(final) #passo da integração por MCMC
    t=t+1
  }
  print("rodando...")
  tabela_resultados[i]=mean(media)
  i=i+1
}


round(tabela_resultados, 2) #resultado final arredondado


