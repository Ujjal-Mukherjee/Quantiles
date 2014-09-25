def matrix(m=1,n=1):
	#Create a zero matrix
	new_matrix = [[0 for row in range(n)] for col in range(m)]
	return new_matrix

x=matrix(10,10)

print x

def dim(x):
	#Dimension of a matrix
	dim = [len(x), len(x[1])]
	return dim

d=dim(x)
print d

def display_matrix(x):
	n=dim(x)[0]
	k=dim(x)[1]
	for i in range(n):
		print x[i]
		


print display_matrix(x)

def matrix_mult(x,y):
	n1=dim(x)[0]
	n2=dim(y)[0]
	k1=dim(x)[1]
	k2=dim(y)[1]

	if(k1!=n2):
		print "non conformable matrix"
		return

	mult_matrix=matrix(n1,k2)

	for i in range(n1):
		for j in range(k2):
			mult=0
			for k in range(k1):
				mult=mult+x[i][k]*y[k][j]
			
			mult_matrix[i][j]=mult

	return mult_matrix

x=matrix(10,20)
y=matrix(20,30)

z=matrix_mult(x,y)

display_matrix(z)

f=open('TestDataPerl.csv', 'r')

x=f.readlines()

f.close()


def read_data(x, sep=','):
	n=dim(x)[0]
	k1=dim(x)[1]

	k=1

	for i in range(k1):
		if(x[1][i]==sep):
			k=k+1
	
	data=matrix(n,k)

	for i in range(n):
		k1=len(x[i])
		counter=0
		temp=''

		for j in range(k1):
			
			if((x[i][j]==sep) or (x[i][j]=='\n') or (x[i][j]=='\r')):
				
				if(counter<k):
					data[i][counter]=temp
					temp=''				
					counter=counter+1
			else:
				
				temp=temp + x[i][j]

	

	return data

data=read_data(x)
print "******************************"
#display_matrix(data)
#print data[0][0]
#print dim(data)
#print "******************************"


def col_vec(data,n,header='T',num='T'):

	n1=dim(data)[0]
	k1=dim(data)[1]

	if(n>k1):
		print 'Not a valid col index'
		return

	if(header=='T'):
		R=range(1,n1)
		R1=range(n1-1)
	else:
		R=range(n1)
		R1=range(n1)

	col_n=[0 for row in R1]

	for i in R:
		if(header=='T'):
			if(num=='T'):
				col_n[i-1]=float(data[i][n])
			else:
				col_n[i-1]=data[i][n]
		else:
			if(num=='T'):
				col_n[i]=float(data[i][n])
			else:
				col_n[i]=data[i][n]

	return col_n


#y=col_vec(data,5)


y=col_vec(data,2,'T','F')

print y

print len(y)

z=[float(i)*float(i) for i in y]

print z

