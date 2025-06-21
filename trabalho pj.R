# Configurações iniciais 

library(parallel)
library(ggplot2)
library(DescTools)
library(gridExtra)
library(scales)
set.seed(333)

# Distrinuições a Serem Testadas (MODIFICAR PARA NORMAL E N DISTRIBUIÇÕES T COM DIFERENTES GRAUAS DE LIBERDADE)

x11("Distribuições Teóricas")
p1_1 <- ggplot(data.frame(x=c(0,1)), aes(x=x)) + 
    stat_function(fun=dnorm, colour="black", size=1.5) +
    xlim(-4,4) + 
    ggtitle("Distribuição: Normal Padrão")

p1_2 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
  stat_function(fun = dt, args = list(df = 50), color = "black", size = 1.5) +
  ggtitle("Distribuição: t (50 graus de liberdade)")

p1_3 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
  stat_function(fun = dt, args = list(df = 40), color = "black", size = 1.5) +
  ggtitle("Distribuição: t (40 graus de liberdade)")

p1_4 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
  stat_function(fun = dt, args = list(df = 30), color = "black", size = 1.5) +
  ggtitle("Distribuição: t (30 graus de liberdade)")

p1_5 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
  stat_function(fun = dt, args = list(df = 20), color = "black", size = 1.5) +
  ggtitle("Distribuição: t (20 graus de liberdade)")

p1_6 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
  stat_function(fun = dt, args = list(df = 10), color = "black", size = 1.5) +
  ggtitle("Distribuição: t (10 graus de liberdade)")

grid.arrange(p1_1, p1_2, p1_3, p1_4, p1_5, p1_6, ncol=3)

# Tamanhos Amostrais
#N <- c(seq(10, 100, by = 10), seq(150,500, by = 50), seq(600,1000, by = 100), seq(2000, 10000, by = 1000))
N <- c(10,100,1000)
# Funções

# 1) Testes não paramétricos
testeNaoParametricos <- function(N){
    
    # MODIFICAR AS DISTRIBUIÇÕES PARA NORMAL E T
    # Amostra Normal
    amostraNormal <- rnorm(N)
  
    # Amostra T, 100 graus de liberdade
    amostraT100 <- rt(N, 100)

    # Amostra T, 80 graus de liberdade
    amostraT80 <- rt(N, 80)

    # Amostra T, 60 graus de liberdade
    amostraT60 <- rt(N, 60)

    # Amostra T, 40 graus de liberdade
    amostraT40 <- rt(N, 40)

    # Amostra T, 20 graus de liberdade
    amostraT20 <- rt(N, 20)
  
    # Kolmogorov-smirnov Test
    kolmogorovNormal <- ks.test(amostraNormal,
                              "pnorm",
                              mean=mean(amostraNormal),
                              sd=sd(amostraNormal))

    kolmogorovT100 <- ks.test(amostraT100,
                              "pnorm",
                              mean=mean(amostraT100),
                              sd=sd(amostraT100))
    kolmogorovT80 <- ks.test(amostraT80,
                              "pnorm",
                              mean=mean(amostraT80),
                              sd=sd(amostraT80))
    kolmogorovT60 <- ks.test(amostraT60,
                              "pnorm",
                              mean=mean(amostraT60),
                              sd=sd(amostraT60))
    kolmogorovT40 <- ks.test(amostraT40,
                              "pnorm",
                              mean=mean(amostraT40),
                              sd=sd(amostraT40))
    kolmogorovT20 <- ks.test(amostraT20,
                              "pnorm",
                              mean=mean(amostraT20),
                              sd=sd(amostraT20))

    kolmogorovResultado <- c(kolmogorovTipoI            = kolmogorovNormal$p.value < 0.05,
                             kolmogorovTipoIIT100       = !(kolmogorovT100$p.value < 0.05),
                             kolmogorovTipoIIT80        = !(kolmogorovT80$p.value < 0.05),
                             kolmogorovTipoIIT60        = !(kolmogorovT60$p.value < 0.05),
                             kolmogorovTipoIIT40        = !(kolmogorovT40$p.value < 0.05),
                             kolmogorovTipoIIT20        = !(kolmogorovT20$p.value < 0.05))

    # Robust Jarque-Bera Test
    jarqueNormal <- JarqueBeraTest(amostraNormal)
    jarqueT100 <- JarqueBeraTest(amostraT100)
    jarqueT80 <- JarqueBeraTest(amostraT80)
    jarqueT60 <- JarqueBeraTest(amostraT60)
    jarqueT40 <- JarqueBeraTest(amostraT40)
    jarqueT20 <- JarqueBeraTest(amostraT20)
  
    jarqueResultado <- c(jarqueTipoI            = jarqueNormal$p.value < 0.05,
                         jarqueTipoIIT100       = !(jarqueT100$p.value < 0.05),
                         jarqueTipoIIT80        = !(jarqueT80$p.value < 0.05),
                         jarqueTipoIIT60        = !(jarqueT60$p.value < 0.05),
                         jarqueTipoIIT40        = !(jarqueT40$p.value < 0.05),
                         jarqueTipoIIT20        = !(jarqueT20$p.value < 0.05))

    # Resultados
    resultado <- c(kolmogorovResultado, jarqueResultado)
    return(resultado)
}

# 2) Replicação Paralela
replicacaoParalela <- function(N, cl) {

    t1 <- Sys.time()

    # Equivalente ao replicate, para paralelismo
    simulacao <- parSapply(cl, 1:1E5, function(i) testeNaoParametricos(N))

    # Calculo da proporção de TRUEs, que representam os erros.
    resultado <- rowMeans(simulacao)

    print(paste0("N = ", N))
    t2 <- Sys.time()
    dif <- t1-t2
    print(dif)

    return(c(N = N, resultado))
}

# Criação do cluster e exportação das funções
nCluSter <- 12 # IMPORTANTE: Substituir pelo número de cores que deseja alocar
cl <- makeCluster(nCluSter)
clusterExport(cl, c("testeNaoParametricos", "N"), envir = environment())
clusterEvalQ(cl, library(DescTools))


# Execução da Replicação e calculo do tempo total
t1 <-  Sys.time() 
resultado <- sapply(N, function(x) replicacaoParalela(x, cl))
t2 <- Sys.time()
t2 - t1

# Encerramento do cluster
stopCluster(cl)

resultado <- t(resultado)

save(resultado, file = "resultado.Rdata")

#load("resultado.Rdata")