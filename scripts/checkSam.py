# Date: Jun 29, 2015 $
# Author: Derek Guan $
# Purpose: check Sam files of simulated data #
# params: 1st: maf 2nd sam files or just sam files 3rd limit 
# params: 1st sam file
import sys,re

def exist(samlist):
	samlen = len(samlist)
	i = 11
	indsa = -1
	indxa = -1
	while i < samlen:
		if samlist[i][0:2] == "SA":
			indsa = i
		elif samlist[i][0:2] == "XA":
			indxa = i
		i = i + 1
	return (indsa,indxa)
def dealPSam(sams, length):
	# sams list of sam
	# includes isrc, cigar
	# return least S, 
	i = 0
	posSign = []
	while i < length:
		posSign.append(0)
		i = i + 1
	min_clip = 0xFFFFFFFF
	for sam in sams:
		# find rclip and lcip
		lclip = 0
		rclip = 0
		m = re.match("\d+[SH]",sam[3])
		if m:
			lclip = int(m.group(0)[:-1])

		m = re.search("\d+[SH]$",sam[3])
		if m:
			rclip = int(m.group(0)[:-1])
		if rclip + lclip < min_clip:
			min_clip = rclip + lclip
		if sam[1] == 16:
			#exchange lclip 
			temp = lclip
			lclip = rclip
			rclip = temp
		i = lclip
		end = length - rclip
		while i < end:
			posSign[i] = 1
			i = i + 1
	coverage = 0
	for e in posSign:
		if e == 1:
			coverage = coverage + 1
	return (coverage,min_clip)	

def dealSam(std, sams):
	#decription of these parameters:
	# std: a list includes several element 
	# 1st Element: chrIndex 
	# 2nd Element: isRC:i:0|16
	# 3rd Element: lenRead
	# 4th Element: pos: position on reference that corresponds to a specific character in read.
	# sams: a list includes several element, each element is a list
	# Element: chrIndex, ISRC, pos, cigar 
	posSign = []

	quaSign = []

	rdlength = std[2]

	posList = std[3]

	i = 0
	while i < rdlength:
		posSign.append(0)
		i = i + 1
	limit = 50
	
	min_diff = 0xFFFFFFFF
	qualified = False
	for sam in sams:
		#input()
		lclip = 0
		rclip = 0
		m = re.match("\d+[SH]",sam[3])
		if m:
			lclip = int(m.group(0)[:-1])

		m = re.search("\d+[SH]$",sam[3])
		if m:
			rclip = int(m.group(0)[:-1])
		if lclip >= rdlength or rclip >= rdlength or rclip + lclip >= rdlength:
			continue
		

		if std[0] == sam[0] and std[1] == sam[1]:
			# find rclip and lcip
			diff = abs(posList[lclip] - sam[2])
			leftLength = rdlength - lclip - rclip

			if diff <= limit: 
				start = lclip
				end = rdlength - rclip
				while start < end:
					mappedSign[start] = 1
					start = start + 1
				if float(leftLength) >= 0.8 * float(rdlength):
					qualified = True

			if diff < min_diff:
				min_diff = diff
			
		# count mapped length, reflect to positive side	

		if sam[1] == 16:
			temp = rclip
			rclip = lclip
			lclip = temp

		start = lclip
		end = rdlength - rclip
		while start < end:
			posSign[start] = 1
			start = start + 1
			
	mappedlen = 0
	for ele in posSign:
		if ele == 1:
			mappedlen = mappedlen + 1
	qualen = 0
	for ele in quaSign:
		if ele == 1:
			qualen = qualen + 1		
	if qualified:
		return ("Y", mappedlen, qualen, min_diff)
	else:
		return ("N", mappedlen, qualen, min_diff)

def AddSams(samlnlist, sams):
	sam_isrc = 0x10&int(samlnlist[1])
	sam_chrname = samlnlist[2]
	sam_pos = int(samlnlist[3])
	sam_cigar = samlnlist[5]
	sams.append([sam_chrname,sam_isrc,sam_pos,sam_cigar])
	(indsa,indxa) = exist(samlnlist)
	#get SA
	#print(indsa)
	if indsa != -1:
		SAs = samlnlist[indsa][5:].strip(";").split(";")
		for SA in SAs:
			salist = SA.split(",") 
			#print(salist)
			sam_chrname = salist[0]
			sam_pos = int(salist[1])
			if salist[2] == '-':
				sam_isrc = 16
			else:
				sam_isrc = 0
			sam_cigar = salist[3]
			sams.append([sam_chrname,sam_isrc,sam_pos,sam_cigar])
	if indxa != -1:
		XAs = samlnlist[indxa][5:].strip(";").split(";")
		for XA in XAs:
			salist = XA.split(",")
			#print(salist)
			sam_chrname = salist[0]
			sam_pos = int(salist[1][1:])
			if salist[1][0] == '-':
				sam_isrc = 16
			else:
				sam_isrc = 0
			sam_cigar = salist[2]
			sams.append([sam_chrname,sam_isrc,sam_pos,sam_cigar])

def main():
	samfl = open(sys.argv[1])
	chrNames = []
	for samline in samfl:
		if samline[0] == '@':
			if samline[1:3] == "SQ":
				chrNames.append(samline.split('\t')[1][3:])
		else:
			samline = samline.strip('\n')
			break
	#chrNames.append("gi|110640213|ref|NC_008253.1|")
	#print(len(sys.argv))
	stdfl = open(sys.argv[2])
	if len(sys.argv) >= 4:
		
		limit = int(sys.argv[3])
		# get all chrosome's Name
		samlnlist = samline.split('\t')
		seq = samlnlist[0].split('/')[0]

		for stdline in stdfl:	
			secstdline = stdfl.readline().strip('\n')
			thrdstdline = stdfl.readline().strip('\n')
			stdfl.readline()

			thrdstdlist = thrdstdline.split()
			secstdlist = secstdline.split()
			rdlength = int(thrdstdlist[5])
			stdpos = int(secstdlist[2])
			
			
			if seq != thrdstdlist[1]:
				#print(seq)
				#print(thrdstdlist[1])
				print("%s\t%s\t%d\t%d\t%d\t%d" %(thrdstdlist[1],"N", rdlength, -1, -1, -1))
				continue
			if samlnlist[5] == "*":
				print("%s\t%s\t%d\t%d\t%d\t%d" %(thrdstdlist[1],"N", rdlength, -1, -1, -1))
				samlnlist = samfl.readline().strip('\n').split('\t')
				seq = samlnlist[0].split('/')[0]
				continue

			strlen = len(thrdstdlist[6])
			q = 0
			poslist = []
			pos = stdpos
			while q < strlen:
				if secstdlist[6][q] == '-':
					poslist.append(pos)
				elif thrdstdlist[6][q] == '-':
					pos += 1
				else:
					poslist.append(pos)
					pos += 1
				q += 1
			#build std
			if thrdstdlist[4] == '-':
				stdsam_isrc = 16
			else:
				stdsam_isrc = 0
			#print(thrdstdlist[1])
			Chrindex = int(thrdstdlist[1].split('_')[0][1:])
			
			stdsam = [chrNames[Chrindex - 1],stdsam_isrc,rdlength,poslist]
			sams = []
			while seq == thrdstdlist[1]:
				#build position list
				# intiate the first 
				AddSams(samlnlist,sams)
				# get next check seq
				samlnlist = samfl.readline().strip('\n').split('\t')
				seq = samlnlist[0].split('/')[0]
				# if the same create next otherwise deal sam 

			(yes,mappedlen, qualen,min_diff)= dealSam(stdsam,sams)	
			print("%s\t%s\t%d\t%d\t%d\t%d" %(thrdstdlist[1],yes, rdlength, mappedlen, qualen, min_diff))
	else:

		samlnlist = samline.split('\t')
		seq = samlnlist[0]

		# get next check seq
		for stdline in stdfl:	
			secline = stdfl.readline().strip('\n')
			stdfl.readline()
			stdfl.readline()
			stdseq = stdline.strip('\n').split()[0][1:]

			rdlength = len(secline)
			# print(stdseq)
			# print(seq)
			if stdseq != seq: 
				print("%s\t%s\t%d\t%d\t%d" %(stdseq,"N", rdlength, -1, -1))
				continue
			if samlnlist[5] == '*':
				print("%s\t%s\t%d\t%d\t%d" %(stdseq,"N", rdlength, -1, -1))
				samlnlist = samfl.readline().strip('\n').split('\t')
				seq = samlnlist[0]
				continue
			sams = []
			while seq == stdseq:
				AddSams(samlnlist,sams)
				samlnlist = samfl.readline().strip('\n').split('\t')
				seq = samlnlist[0]
					#rdlength = len(newlist[9])
			(mappedlen, min_clip) = dealPSam(sams,rdlength)
			#print(sams)
			print("%s\t%s\t%d\t%d\t%d" %(stdseq, "Y", rdlength, mappedlen, min_clip))

if __name__ == '__main__':
	main()