

# observation 
obs = read.table('../data/obs.csv', header=T, sep=',')
plot(obs, pch = 19, ylab = '# of New infections')

## simulation 
## 5d input
d <- matrix(scan('data/sim_design.txt'), ncol=5, byrow=T)
colnames(d) <- paste0('theta_', 1:5)

pairs(d, gap=.01,pch=16,cex=.5, col = 'grey',
      labels=c(expression(theta[1]), expression(theta[2]), expression(theta[3]), 
                                                expression(theta[4]), expression(theta[5])))

## time-series output
load('data/sim.RData')
nrep = 100; 
nweek = 56; 
m = 100

## sima is an nrep x m x nweek dimensional array consisting of 
## all simulations. The simulation corresponding to
## i-th row in d (i.e. d[i,]) can be obtained by 
## calling sima[,i,].

## Plot the simulations
par(mfrow=c(10,10),oma=c(4.1,5.1,1,1),mar=c(0,0,0,0))
for(k in 1:100){
  matplot(1:nweek, t(sima[,k,])+1, type='l', ylim=range(sima+1), axes=F, log='y', col="grey")
  box() 
  if(k%%20==1) axis(2, las = 1)
  if(k > 90 & k%%2==1) axis(1)
}
mtext('cumulative number infected',side=2,line=4,outer=T)
mtext('week',side=1,line=3,outer=T)

