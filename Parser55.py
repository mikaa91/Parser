# encoding utf 8 
import numpy
import string


class parser(object):
	def __init__(self,typ,plikin=None): 
		self.npop=0
		self.nloc=0
		self.osinpop=[]
		self.allele=[] # allele==[  pop1[  os1[ l1(A,a),  l1(B,B),  l3(c,c)  ],  os2[]  ],  pop2[],  pop3[]  ]
		self.popids=[]
		self.locids=[]
		self.alinloc=[]
		self.descrip=''
		self.maxallsig=0
		self.datarepr=1
		if plikin==None: plikin=input('Input file name?\n')
		tresc=readfile(plikin)
		if typ=="arlequin":
			self.descrip=tresc[1].split('"')[0]
			i=2
			while i<len(tresc):
				if "SampleName" in tresc[i]:
					self.npop+=1
					self.popids.append(tresc[i].strip().split('"')[1])
				elif "SampleSize" in tresc[i]:
					self.osinpop.append(int(tresc[i].strip().split('=')[1]))
				elif "SampleData" in tresc[i]:
					i+=1
					if self.nloc==0: self.nloc=len(tresc[i].split())-2
					allele=[]
					for j in range(self.osinpop[self.npop]):
						allele.append([i] for i in tresc[i].split()[2:])
						i+=1
						for k in range(self.nloc):
							allele[k].append(tresc[i].split()[k])
						self.allele[self.npop].append([tuple(i) for i in allele])
						i+=1
				i+=1

		elif typ=="linkdos":
			self.npop=int(tresc[0].split()[0])
			self.nloc=int(tresc[0].split()[1])
			self.allele=[[] for i in range(self.npop)]
			i=0
			for j in range(self.npop*2):
				i+=1
				if i%2==1: self.popids.append(tresc[i]) #popid
				else: self.osinpop.append(int(tresc[i]))# osinpop
			while len(tresc[i])==len(tresc[i].lstrip()):
				self.locids.append(tresc[i].split()[0])
				self.alinloc.append(tresc[i].split()[1])
				i+=1
			for j in range(len(self.osinpop)): #pop
				for k in range(self.osinpop[j]): #osob
					linia=tresc[i].split()
					self.allele[j].append([(linia[i],linia[i+1]) for i in range(0,len(linia),2)])
					i+=1
			self.maxallsig=max(self.alinloc)



		elif typ=="fstat":
			self.npop=int(tresc[0].split()[0])
			self.nloc=int(tresc[0].split()[1])
			self.maxallsig=int(tresc[0].split()[2])
			self.datarepr=int(tresc[0].split()[3])
			self.allele=[[] for i in range(self.npop)]
			i=1
			while len(tresc[i].split())==1:
				self.locids.append(tresc[i])
				i+=1
			while i<len(tresc): #n columns of locis + column of pop nb
				linia=tresc[i].split()
				if len(self.osinpop)<int(linia[0]): self.osinpop.append(0)
				self.osinpop[int(linia[0])-1]+=1
				self.allele[int(linia[0])-1].append([(k[:2],k[2:]) for k in linia[1:]])
				i+=1
			self.popids=[str(i+1) for i in range(self.npop)]
			self.alinloc=[0]*self.nloc 


		elif typ=="genepop":
			self.descrip=str(tresc[0])
			self.locids=[i for i in tresc[1].strip().replace(' ','').split(',')]
			self.nloc=len(self.locids)
			
			i=2 #1st"pop"
			while i<len(tresc):
				if 'pop' in tresc[i].lower():
					self.npop+=1
					self.allele.append([])
					self.osinpop.append(0)
					jestpop=False #get out popid if in
				else:
					if not jestpop: #if there's no those popid
						self.popids.append(tresc[i].split(',')[0])
						jestpop=True
					linia=tresc[i].split()
					self.osinpop[self.npop-1]+=1
					self.allele[self.npop-1].append([(k[:2],k[2:]) for k in linia[1:]])
				i+=1
			self.alinloc=[0]*self.nloc



	def zapisz(self,typ,plikout=None):
		if plikout==None: plikout=input('Output file name?\n')
		wynik=open(plikout,'w')
		if typ=="arlequin":
			header='[Profile]\n\tTitle="%s"\n\n\tNbSamples=%s\n\tDataType=STANDARD\n\tGenotypicData=1\n\tLocusSeparator=WHITESPACE\n\tGameticPhase=0\n\tRecessiveData=0\n\tMissingData="?"\n\n[Data]\n\t[[Samples]]\n' % (self.descrip,self.npop)
			end='}'
			sep=" "
			wynik.writelines(header)
			for i in range(self.npop):
				wynik.writelines('\t  SampleName="%s"\n\t  SampleSize=%s\n\t  SampleData= {\n' % (self.popids[i],self.osinpop[i]))
				for j in range(self.osinpop[i]):
					wynik.writelines(str(j+1)+' '*5+'1'+' ')
					for k in range(self.nloc):
						wynik.writelines(' '*2+str('%.2d' % int((self.allele[i][j][k][0]))))
					wynik.writelines('\n'+' '*8)
					for k in range(self.nloc):
						wynik.writelines(' '*2+str('%.2d' % int((self.allele[i][j][k][1]))))
					wynik.writelines('\n')
				wynik.writelines('\t  '+end+'\n')
			
		elif typ=="fstat":
			sep=' '
			wynik.writelines("%s %s %s %s" % (self.npop,self.nloc,self.maxallsig,self.datarepr)+"\n")
			for i in self.locids: wynik.writelines(i+'\n')
			for i in range(self.npop):
				for j in self.allele[i]:
					a=[]
					for k in j:
						a.extend(k)
					wynik.writelines(self.popids[i]+5*sep)
					for k in range(len(a)-1):
						wynik.writelines(a[k])
						if k%2==1: wynik.writelines(sep)
					wynik.writelines(a[-1]+"\n")
			wynik.close()

		elif typ=="genepop":
			sep1=","
			sep2=" "
			wynik.writelines(self.descrip + "\n     ")
			for i in range(self.nloc-1):
				wynik.writelines(self.locids[i] +sep1 + sep2)
			wynik.writelines(self.locids[-1]+'\n')
			for i in range(self.npop):
				wynik.writelines("POP\n")
				for j in self.allele[i]:
					wynik.writelines(self.popids[i]+sep1)
					a=[]
					for k in j:
						a.extend(k)
					wynik.writelines(sep2)
					for k in range(len(a)):
						wynik.writelines(str('%.2d' % int(a[k])))
						if k%2==1: wynik.writelines(sep2)
					wynik.writelines("\n")

		elif typ=="linkdos":
			sep="\t"
			wynik.writelines("%s %s" % (self.npop,self.nloc)+"\n")
			for i in range(self.npop):
				wynik.writelines(self.popids[i]+"\n"+str(self.osinpop[i])+"\n")
			for i in range(self.nloc):
				wynik.writelines(self.locids[i]+sep+str(self.alinloc[i])+"\n")
			for i in self.allele:
				for j in i:
					a=[]
					for k in j:
						a.extend(k)
					for k in a:
						wynik.writelines(sep+str(int(k)))
					wynik.writelines("\n")
		wynik.close()




def readfile(nazwa):
	plik=open(nazwa)
	inputf=plik.readlines()
	plik.close()
	for i in range(len(inputf)):
		inputf[i]=inputf[i].rstrip('\n').rstrip('\r').rstrip('\n')
	return inputf


typin=input('Choose input data format: \nfstat\ngenepop\nlinkdos\narlequin\n')
typout=input('Choose required output format:\nfstat\ngenepop\nlinkdos\narlequin\n')
parserroboczy=parser(typin)
parserroboczy.zapisz(typout)














