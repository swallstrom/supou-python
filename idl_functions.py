import numpy as np

#percentile
### This can be done with np.percentile(data,alpha)
### Also, not used in the supou code...
def percentile(data,alpha,index=index):

	sorted = np.sort(data)

	ndata = np.shape(data)[0]

	index = alpha*ndata
	index = sorted[index]

	return, data[index]

#permute
### Not used in the supou code....
def permute(n,seed):

	unif = np.randomu(seed,n)

	return, ind[np.sort(unif)]


#logit and inv_logit
def logit(x,inverse=False):

	if inverse == False:
		xvals == x[(x < 0) or (x > 1)]

		if True in xvals: 
			print('All elements of x must be between 0 and 1')
			return, 0

		logit = np.log(x / (1 - x)) #Is this natural log? Yes.

		return, logit

	else:
		inv_logit = exp(x) / (1 + exp(x))

		return, inv_logit

#gettok
### Not used in the supou code...
def gettok(st,char,exact=True,notrim=True):


#character is a blank, treat tabs as blanks
	if exact == True:
		st.strip()
		if char == ' ':
			for i in range(0,len(st)):
				test = st[i].find("\t")

				if test != -1: 
					st[i].replace("\t",char)

				st[i].split(char)

	token = st[i][0]

#find character in string
	for i in range(0,len(st)):
		pos = st[i].find(char)
		if pos == -1: 
			st[i] = ''
			print(i,st[i])
		if pos != -1: 
			st[i] = st[i].split(char)
			st[i] = st[i][0]
			print(i,st[i])
	
#return result
	return, token

#ntostr
### This function is just used to convert a float to string so it can be printed out... use str(num) instead
def ntostr(num,pos2,pos1,format=format,round=False):
	
	strnum = []
	nshape = np.shape(c)[0]
#because it would be too fucking hard for Python to work on arrays
#of strings, right?
	for i in range(0,nshape):
		strnum.append(str(num[i]))

		if len(pos1) == 0: pos1 = 0
		if len(pos2) == 0: pos2 = len(str(num[i]))

#default return whole string unless otherwise specified
		if round:
			nleft = long(alog10(abs(num[i])))+1 > 1
			nleft = nleft + (num[i] < 0)
			nav = pos2 - nleft - 1

			if nav >= 0:
				tnum = rnd(num[i],nav)
			else:
				nav = nav+1
				tnum = rnd(num,nav)
				
			strnum[i] = str(tnum)


		else:
			strnum[i] = strnum[i][pos1:pos2+1]

	return,strnum

#function to return value rounded at required precision
def rnd(value,round,direction):
	fnctn = value*0.0
	nonzero = value[value != 0.]

	if nonzero > 0 :

		if (n_params() lt 2) then round = 1.0
		round = abs(round)
		val=value[nonzero]+round*0.5

		if direction == 'up':
			val = value(nonzero)+round*0.999999 
		elif direction == 'down':
			val = float(long((value[nonzero]+round)/float(round))*round)-round*0.999999
		elif direction == 'out':
			val = value[nonzero]+float(value[nonzero])/abs(value[nonzero])*round*0.999999
		elif direction == 'in':
			val = value[nonzero]-float(value[nonzero])/abs(value[nonzero])*round*0.000001

		fnctn[nonzero] = float(long(val/float(round))*round)

	sz = value.shape()[0]
	if sz[0] == 0:
		fnctn = fnctn[0] #data originally a scalar

	return, fnctn
