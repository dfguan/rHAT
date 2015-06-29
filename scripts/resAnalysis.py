import sys

string=sys.argv[1]
fp = open(sys.argv[1])
lim = int (sys.argv[2])

fp.readline()
strlist = string.split('_')

readcount = 0
error = 0
tlen = 0
erlen = 0
notfound = 0
for eachline in fp:
	readcount += 1
	eachlist = eachline.strip('\n').split('\t')
	diff 	= int(eachlist[7])
	length 	= int(eachlist[3])
	tlen += length
	if diff > lim :
		error += 1
		erlen += length
	elif diff == -1 :
		notfound += 1 	
	
print ("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d" %(strlist[1],strlist[2],strlist[3],strlist[4],readcount,notfound,error,erlen,tlen))
