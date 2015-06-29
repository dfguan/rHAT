# Author: Derek Guan
# Date : Jan 28, 2015

import sys


def dealP(result):
	#notfound=0
	error = 0
	mapped = 0
	notfound = 0
	genNum = 0
	quaGenNum = 0
	
	for l in result:
		llist = l.strip('\n').split('\t')
		genNum = genNum + int(llist[2])
		if llist[1] == 'N':
			notfound = notfound + 1
		else:
			quaGenNum = quaGenNum + int(llist[3])
			mapped = mapped + 1
			if float(llist[4]) > 0.2*float(llist[2]):
				error = error + 1
	return (mapped, notfound,  error, quaGenNum, genNum)


def dealS(result):
	notfound = 0
	error = 0
	readNum = 0
	genNum = 0
	quaGenNum = 0
	mappedGenNum = 0
	for l in result:
		llist = l.strip('\n').split('\t')
		readNum = readNum + 1
		genNum = genNum + int(llist[2])
		
		if llist[1] == 'N':
			if llist[3] == '-1':
				notfound = notfound + 1
			else:
				error = error + 1
				mappedGenNum = mappedGenNum + int(llist[3])
				quaGenNum = quaGenNum + int(llist[4])

		else:
			mappedGenNum = mappedGenNum + int(llist[3])	
			quaGenNum = quaGenNum + int(llist[4])

	return (readNum, notfound, error, quaGenNum, genNum)

def main():
	resultName = sys.argv[1]
	result = open(resultName)
	resNameList = resultName.strip('\n').split('_')
	lenresNameList = len(resNameList) 
	if lenresNameList == 9:
		# bwa blasr pacbio data
		(mapped, notfound, error, quaGenNum, genNum) = dealP(result)
		print("%s\t%s\t%d\t%d\t%d\t%d\t%d" %(resNameList[1],resNameList[2],mapped,notfound, error,quaGenNum,genNum))
	elif lenresNameList == 5:
		#bwa blasr simulated data
		(readNum, notfound, error, quaGenNum, genNum) = dealS(result)
		print("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d" %(resNameList[1],resNameList[2],resNameList[3],resNameList[4],readNum,notfound,error,quaGenNum,genNum))
	elif lenresNameList == 11:
		# exp result
		#exp pacbio data
		(mapped,notfound, error, quaGenNum, genNum) = dealP(result)
		print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d" %(resNameList[8],resNameList[9],resNameList[10],mapped, notfound, error, quaGenNum,genNum))
	else:
		# exp simulated data
		(readNum, notfound, error, quaGenNum, genNum) = dealS(result)
		print("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d" %(resNameList[1],resNameList[2],resNameList[3],resNameList[4],resNameList[5],readNum,notfound,error,quaGenNum,genNum))
	
if __name__ == '__main__':
	main()	