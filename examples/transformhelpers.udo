
#define DEFAULT_SIZE #1024#
#define DEFAULT_HOPSIZE #128#

opcode tfin, k[], ajj
	a1, isize, ihopsize xin
	if (isize == -1) then
		isize = $DEFAULT_SIZE
	endif
	if (ihopsize == -1) then
		ihopsize = $DEFAULT_HOPSIZE
	endif
	
	iolaps = isize/ihopsize
	kcnt init 0
	krow init 0
	kTemp[] init isize
	kOla[] init isize
	kIn[] init isize
	kOut[][] init iolaps, isize 

	if (kcnt == ihopsize) then  
		kWin[] window kIn, krow*ihopsize
		kSpec[] tfhaar1 kWin
		krow = (krow+1) % iolaps
		kcnt = 0
	endif

	kIn shiftin a1

	xout kSpec
	kcnt += ksmps
endop


opcode tfout, a, k[]jj
	kSpec[], isize, ihopsize xin
	if (isize == -1) then
		isize = $DEFAULT_SIZE
	endif
	if (ihopsize == -1) then
		ihopsize = $DEFAULT_HOPSIZE
	endif

	iolaps = isize/ihopsize ; overlaps
	kcnt init 0    ; counting vars
	krow init 0
	kTemp[] init isize
	kOla[] init isize ; overlap-add buffer
	kIn[] init isize  ; input buffer
	kOut[][] init iolaps, isize ; output buffers

	if (kcnt == ihopsize) then  
		kRow[] tfhaar1inv kSpec
		kWin[] window kRow, krow*ihopsize
		kOut setrow kWin, krow
		kOla = 0
		ki = 0
		until (ki == iolaps) do
			kRow getrow kOut, ki
			kOla = kOla + kRow
			ki += 1
		od
 
		krow = (krow+1)%iolaps
		kcnt = 0
	endif
	a2 shiftout kOla
	aout = (a2/iolaps)
	xout aout 
	kcnt += ksmps
endop


opcode scramble, k[], k[]
	kin[] xin
	inum random lenarray(kin)*0.2, lenarray(kin)
	kndx = 0
	while (kndx < inum) do
		ksrc random 0, lenarray(kin)
		kdest random 0, lenarray(kin)
		kin[int(kdest)] = kin[int(ksrc)]
		kndx += 1
	od
	xout kin
endop


opcode invert, k[], k[]
	kin[] xin
	kout[] init lenarray(kin)
	kwritex = 0
	imax = lenarray(kin) - 1
	while (kwritex < lenarray(kin)) do
		kout[kwritex] = kin[imax - kwritex]
		kwritex += 1
	od
	xout kout
endop
