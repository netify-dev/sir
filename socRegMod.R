####
source(paste0(here::here(), '/setup.R'))
####

####
# Load data
load(paste0(inPath, 'socRegData.rda'))

# Create binary joint democracy measure
iPol = ifelse(Z[,,'i_polity',]>17,1,0)
jPol = ifelse(Z[,,'j_polity',]>17,1,0)
jointDem = iPol * jPol
Z[,,'i_polity',] = jointDem

# spec covars
wVars = c(
	'int', 'exports', 'ally',
	'verbCoop','minDistLog'
	)
W = W[,,wVars,] ; W_avg = W_avg[,,wVars]
W[,,c('verbCoop'),] = log(W[,,c('verbCoop'),] + 1)
W_avg[,,c('verbCoop')] = log(W_avg[,,c('verbCoop')] + 1)
zVars = c(
	'int' , 'mConf', 'mConf_ji',
	'minDistLog', 'i_polity',
	'ally', 'exports',
	'verbCoop', 'verbCoop_ji'
	)
Z = Z[,,zVars,]
####

####
# run mods

# run als approach with varying W
start = Sys.time()
dynMod = poisblr_dyn(Y, W, X, Z)
end = Sys.time()
print(end-start)
save(start, end, file=paste0(outPath, 'alsModResults_time.rda'))

# run als approach with static W
staticMod = poisblr(Y, W_avg, X, Z)

# Compare with simple dyadic model
dyadMod = summary(glm( c(Y) ~ -1 + apply(Z,3,c) ,family="poisson" ))$'coefficients'
rownames(dyadMod) = gsub('apply(Z, 3, c)', '(Z) ', rownames(dyadMod), fixed=TRUE)
predGLM = exp(t(dyadMod[,1] %*% t(apply(Z,3,c))))
dyadMod = list(summ=dyadMod, pred=predGLM)
####

####
# clean up var names
varKey=data.frame(dirty=rownames(dynMod$summ),clean=rownames(dynMod$summ),
	stringsAsFactors = FALSE)
for(v in paste0('(',c('Z','betaW','alphaW'),') ') ){
	varKey$clean = gsub(v,'',varKey$clean, fixed=TRUE) }
varKey$cleanTeX = varKey$clean
varKey$cleanTeX[varKey$cleanTeX=='int'] = 'Intercept'
varKey$cleanTeX[varKey$cleanTeX=='mConf'] = 'Material\nConflict$_{ij,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='mConf_ji'] = 'Material\nConflict$_{ji,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='minDistLog'] = 'Distance$_{ij,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='verbCoop'] = 'Verbal\nCooperation$_{ij,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='verbCoop_ji'] = 'Verbal\nCooperation$_{ji,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='i_polity'] = 'Joint\nDemocracy$_{ij,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='exports'] = 'Log(Trade)$_{ij,t-1}$'
varKey$cleanTeX[varKey$cleanTeX=='ally'] = 'Ally$_{ij,t-1}$'

# org coefData
coefData = list(socInflDyn=dynMod$summ,
	socInflStatic=staticMod$summ, glm=dyadMod$summ)
coefData = lapply(coefData, function(x){
	x = data.frame(x, stringsAsFactors = FALSE)

	# clean up var names
	x$var = rownames(x) ; rownames(x) = NULL
	x$varClean = varKey$cleanTeX[match(x$var, varKey$dirty)]

	# define types
	x$type = 'Z'
	x$type[grepl('(betaW) ',x$var,fixed=TRUE)] = 'betaW'
	x$type[grepl('(alphaW) ',x$var,fixed=TRUE)] = 'alphaW'

	# remove intercept
	x = x[x$varClean !='Intercept',]

	# order var
	x$varClean = factor(x$varClean,
		levels=rev(unique(intersect(varKey$cleanTeX,x$varClean))))

	return(x) })
####

####
# lazy save
save(
	dynMod, staticMod, dyadMod, varKey, coefData,
	file=paste0(outPath, 'alsModResults.rda')
	)
####
