#!MC 1120
# Created by Tecplot 360 build 11.2-0-563
#
#

$!VarSet |MFBD| = '/home/andreas/Codes/particles/gitvringD/VR2D'

$!VARSET |NIT| = 1
$!VARSET |NID| = 1
$!VARSET |NCS| = 5

$!VARSET |IIT|=4000
$!VARSET |IID|=1


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	12
	0.5
	0.6
	0.7
	0.8
	0.9
	1.0
	1.1
	1.2
	1.3
	1.4
        1.5
        20.0
        $!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = BANDED}


	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=2


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat"'

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.5
	0.6
	0.7
	0.8
	0.9
	1.0
	1.1
	1.2
	1.3
	1.4
        1.5
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.5}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.5}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=3


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.5
	0.6
	0.7
	0.8
	0.9
	1.0
	1.1
	1.2
	1.3
	1.4
        1.5
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.5}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.5}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=4


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.0
	0.2
	0.4
	0.6
	0.8
	1.0
	1.2
	1.4
	1.6
	1.8
        2.0
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 2.0}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=5


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.0
	0.2
	0.4
	0.6
	0.8
	1.0
	1.2
	1.4
	1.6
	1.8
        2.0
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 2.0}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=6


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.0
	1.0
	2.0
	3.0
	4.0
	5.0
	6.0
	7.0
	8.0
	9.0
        10.0
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX =10.0}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!VARSET |IIT|=4000
$!VARSET |IID|=7


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.0
	1.0
	2.0
	3.0
	4.0
	5.0
	6.0
	7.0
	8.0
	9.0
        10.0
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX =10.0}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP

$!VARSET |IIT|=4000
$!VARSET |IID|=8


$!READDATASET  '"|MFBD|/post|IIT%5.5d|size|IID%2.2d|.dat" '

  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12" "V13" "V14" "V15" "V16" "V17" "V18" "V19"'


$!LOOP |NCS|
     $!VARSET |IIC| = |LOOP|




      $!IF |IIC| == 1 
	$!GLOBALCONTOUR 1  VAR = 5
      $!ENDIF 
      $!IF |IIC| == 2 
	$!GLOBALCONTOUR 1  VAR = 7
      $!ENDIF 
      $!IF |IIC| == 3 
	$!GLOBALCONTOUR 1  VAR = 8
      $!ENDIF 
      $!IF |IIC| == 4
	$!GLOBALCONTOUR 1  VAR = 9
      $!ENDIF 
      $!IF |IIC| == 5 
	$!GLOBALCONTOUR 1  VAR = 10
      $!ENDIF 




	$!CONTOURLEVELS RESETTONICE
	  CONTOURGROUP = 1
	  APPROXNUMVALUES = 51
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0499999999999996836}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 1.7999999999999996}}
	$!CONTOURLEVELS NEW
	  CONTOURGROUP = 1
	  RAWDATA
	11
	0.0
	1.0
	2.0
	3.0
	4.0
	5.0
	6.0
	7.0
	8.0
	9.0
        10.0
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = 0.0}}
	$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX =10.0}}

	$!EXPORTSETUP EXPORTFORMAT = JPEG
	$!EXPORTSETUP IMAGEWIDTH = 790
	$!EXPORTSETUP EXPORTFNAME = '/home/andreas/Codes/particles/gitvringD/VR2D/particlest|IIT%5.5d|id|IID%2.2d|c|IIC%2.2d|.jpg'
	$!EXPORT 
	  EXPORTREGION = CURRENTFRAME
$!ENDLOOP
$!RemoveVar |MFBD|
