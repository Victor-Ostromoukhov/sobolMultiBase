(* sobolMultiBase.m 1/1/2020
*)
If[(*$ProcessorCount == 4 &&*) Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2]];

SetDirectory[ToFileName[$HomeDirectory,"sobolMultiBase"]];

SetOptions[Graphics, ImageSize ->{ 650,Automatic},AspectRatio->Automatic, PlotRange->All,Frame->False ];
SetOptions[ListDensityPlot, ImageSize -> {512,Automatic},AspectRatio->Automatic, PlotRange->All];
SetOptions[ListPlot, Joined->True, ImageSize ->{1024,Automatic}, PlotStyle->{Red,Blue,Green,Gray}];
SetOptions[StackedListPlot, Joined->True, ImageSize ->{1024,Automatic}, PlotStyle->{Red,Blue,Gray,Green,Cyan,Magenta}];
SetOptions[ListStepPlot, Joined->True, ImageSize ->{1024,Automatic}, PlotStyle->{Red,Blue,Gray,Green,Cyan,Magenta}];
SetOptions[ListLinePlot, Joined->True, ImageSize ->{1024,Automatic}, PlotStyle->{Red,Blue,Gray,Green,Cyan,Magenta}];
Print["module sobolMultiBase loaded."];

(*Needs["TetGenLink`"];*)

pid := "_pid"<>ToString[$ProcessID]<>"_kid"<>ToString[$KernelID]
systemID = StringSplit[$System][[1]]; (* "Mac" or "Linux" *)
(*execPrefix = "/usr/bin/";*)
execPrefix = "LD_LIBRARY_PATH=/usr/local/lib/ "<>"~/bin/";
execPrefix = "~/bin/";

absoluteTime = Floor[100*AbsoluteTime[]];

first := First
second[x_]:= If[Length[x] > 1, x[[2]], First[x] ] (* like First *)
third[x_]:= If[Length[x] > 2, x[[3]], First[x] ] (* like First *)
fourth[x_]:= If[Length[x] > 3, x[[4]], First[x] ] (* like First *)
last := Last

epsilon = 10^-10.;
eps = 10^-6.;

mf := MatrixForm
T :=  Transpose
PI = Pi//N;
known := ValueQ
i2s[n_,len_:6] := ToString[NumberForm[n, len-1, NumberPadding -> "0"]]
n2PaddedString[n_,len_:6] := ToString[NumberForm[n, len-1, NumberPadding -> "0"]]
f2PaddedString[n_,n1_:3,n2_:2] := ToString[NumberForm[n, {n1,n2}, NumberPadding -> "0"]]

i2padded[n_,len_:6] := ToString[NumberForm[n, len-1, NumberPadding -> " "]]

tab2s[tab_] := StringJoin @ Table[" "<>ToString[tab[[i]]]<>" ",{i,Length[tab]}]

euclidlen[z_] := Sqrt[Total[z^2]]
euclidlen2[z_] := Total[z^2]
euclidlenSq[z_] := Total[z^2]
euclidlenN[z_] := Sqrt[Total[z^2]]//N

order2permut[s_] := (Flatten[Drop[#,1]& /@ Sort[Table[{s[[i]],i},{i,Length[s]}]]]); (* 1..n *)
permut2order := order2permut

order2permut1toN := order2permut
order2permut0toNminus1[s_] := ((s+1)//order2permut)-1

colTable = {Red,Green,Blue,Cyan,Magenta,Yellow,LightGray,Gray};
getColor[ind_]:=colTable[[Mod[ind,Length[colTable],1] ]]

niceRaster[img_,OptionsPattern[]] :=
	Block[ {sx,sy,z},
		z = OptionValue[zoom];
		{sy,sx} = Dimensions[img][[;;2]];
		Return[Graphics[Raster[img],PlotRange->{{0,sx},{0,sy}},ImageSize->{z sx,z sy}]];
	];
Options[niceRaster] = {zoom->1};

(*------------------------- supprot for stdSobol -------------------------*)
a058947 = Get["data/a058947.dat"];
(*sobolDirectionVectors = Get["data/sobol_direction_vectors.dat"]; (* 1111 direction vectors taken from new-joe-kuo-6.21201.txt http://web.maths.unsw.edu.au/~fkuo/sobol/
							VO: extra first entry is added to match Kuo,s implementation (by default, parameters direction vectors 1 and 2) *)
*)
sobolDataJoeKuo1111 = Drop[#,1]& @ Import["data/old-joe-kuo-1111.dat"];
sobolDataJoeKuo21201 = Drop[#,1]& @ Import["data/new-joe-kuo-6.21201.dat"];
sobolDataVO = Drop[#,1]& @ Import["data/sobol_init_tab_vo.dat"];

sobolDirectionVectorsJoeKuo = Drop[#,3]& /@ sobolDataJoeKuo21201;

sobolData = Drop[#,1]& @ Import["data/new-joe-kuo-6.21201.dat"];
sobolDirectionVectors = Drop[#,3]& /@ sobolData;

getSobolDegree[sobolInd_]:= Length[IntegerDigits[a058947[[sobolInd]]]]

sobolseqLength = 48; (* to be compatibe with Keller & Grunschloss http://gruenschloss.org/ *)
sobolseqLength = 64;
sobolseqLength = 14;
sobolseqLength = 32;

readKoopmanData[degree_:4] :=
    Module[ {},
    	fname = "data/Koopman/"<>ToString[degree]<>".dat";
    	Print["Reading ",fname];
     	file = OpenRead[fname]; 
   		If[degree <= 32,
    		Do[str = Read[file, String],{4}];
    		data = ToString /@ ReadList[file, String];
	    	len = 1000;
	        numbers = Interpreter["HexInteger"] /@ (RandomSample[data, len]);
	        digits = Join[IntegerDigits[#,2,degree],{1}]& /@ numbers;
	        Print[mf @ digits ];
	        outfname = "data/primitivePolynomialsDegree"<>ToString[degree]<>"_digits_sel"<>ToString[len]<>".dat";
	        Print["Writing into ",outfname];
	        Export[outfname,digits];
	    ,(*ELSE*)
    		data = ToString /@ ReadList[file, String];
	    	len = Length[data];
	        numbers = Interpreter["HexInteger"] /@ (RandomSample[data, len]);
	        numbers = Interpreter["HexInteger"] /@ data;
	        digits = Join[IntegerDigits[#,2,degree],{1}]& /@ numbers;
	        Print[mf @ digits ];
	        outfname = "data/primitivePolynomialsDegree"<>ToString[degree]<>"_digits_sel"<>ToString[len]<>".dat";
	        Print["Writing into ",outfname];
	        Export[outfname,digits];
    	];
    ] (* readKoopmanData *)
    
grayCode[i_, n_] := 
 FromDigits[
  BitXor @@@ Partition[Prepend[IntegerDigits[i, 2, n], 0], 2, 1], 2]
fromGrayCode[i_, n_] := FromDigits[BitXor[IntegerDigits[i, 2, n], FoldList[BitXor, 0, Most[IntegerDigits[i, 2, n]]]], 2]

grayCodeList[k_] := Module[{b = IntegerDigits[k, 2], i},
  Do[
   If[b[[i - 1]] == 1, b[[i]] = 1 - b[[i]]], {i, Length[b], 2, -1} ];
	b
  ]

sobol1DGrayCode[nf_,n_] := (* + Gray Code as described in BratleyFox88 and AntonovSaleev79 *)
    Module[ {i,seq = grayCodeList[n] //Reverse,seqlen,nn},
        seqlen = Length[seq];
        nn = BitXor @@ Table[ seq[[i]] msobol[[nf,i]] 2^(seqlen-i), {i, seqlen}];
        Return[ FromDigits[IntegerDigits[nn,2,seqlen] ,2] / 2^seqlen ]
    ]

buildMSobol[indtab_,dbg_:True,modifiedDirectionVectors_:{},highlight_:0] :=
    Module[ {sobolDataFname,msobol,sobolMX},
    	If[!known[sobolData], 
    		sobolDataFname = "data/sobol_init_tab.dat";  (* a copy of "data/new-joe-kuo-6.21201.dat";*)
        	If[!FileExistsQ[sobolDataFname],Print[sobolDataFname," does not exist"];Abort[] ];
        	If[dbg, Print["loading sobolData from ",sobolDataFname] ];
        	sobolData = Drop[#,1]& @ Import[sobolDataFname];
        ];
        If[ modifiedDirectionVectors =!= {} && Length[indtab] != Length[modifiedDirectionVectors],
            Print["Length[indtab] != Length[modifiedDirectionVectors]"]; Abort[]
        ];
        {msobol,sobolMX} = Table[
        	If[indtab[[i]] == 0,
        		{Table[1,{sobolseqLength}], IdentityMatrix[sobolseqLength]}
        	,(*ELSE*)       	
	            If[ modifiedDirectionVectors === {},
	                getMsobol1D[ indtab[[i]], dbg,{},highlight],(*ELSE*)
	                getMsobol1D[ indtab[[i]],dbg,modifiedDirectionVectors[[i]],highlight ]
	            ]
        	]
        ,{i,Length[indtab]}]//T;
        (*If[ dbg, Print["buildMSobol:  msobol=",msobol//mf]; ];*)
        {msobol,sobolMX}
    ]

getMsobol1D[ind_:0,dbg_:False,modifiedDirectionVectors_:{},highlight_:-1] :=
    Module[ {j,msobol,thisSetDirectionVectors,sobolMX,polynomialDegree,polynomial},
        msobol = Table[1,{sobolseqLength}];
        polynomialDegree = sobolData[[ind,2]];
        polynomial = Join[{1}, IntegerDigits[#,2,polynomialDegree-1]& @ sobolData[[ind,3]], {1} ];

        If[modifiedDirectionVectors === {},
			thisSetDirectionVectors = sobolData[[ind,4;;]]; (* initialization of direction vectors as in new-joe-kuo-6.21201.txt http://web.maths.unsw.edu.au/~fkuo/sobol/ *)
        ,(*ELSE*)
			thisSetDirectionVectors = modifiedDirectionVectors;
        ];
        msobol[[;;polynomialDegree]] = thisSetDirectionVectors;
        Do[msobol[[i]] = BitOr[msobol[[i]],1],{i,sobolseqLength}];
        Do[
            msobol[[i]] = BitXor @@ Join[Table[2^(j) polynomial[[j+1]] msobol[[i-j]],{j,1,polynomialDegree}],{msobol[[i-polynomialDegree]]}];
        ,{i,polynomialDegree+1, sobolseqLength }];
        sobolMX = Table[0,{sobolseqLength},{sobolseqLength}];
        Do[ sobolMX[[i,;;i]] = (IntegerDigits[#,2,i]& @ msobol[[i]])  ,{i, sobolseqLength }];
        If[dbg, getsobolMxGL[msobol,sobolMX,ind,thisSetDirectionVectors,highlight]//Print];
        Return[{msobol,sobolMX}]
    ] (* getMsobol1D *)

getsobolMxGL[msobol_,sobolMX_,ind_,thisSetDirectionVectors_,highlight_:-1,sz_:400] :=
    Module[ {ticks,mm0,mm1,mmm,lbl,polynomialDegree},
    	polynomialDegree = sobolData[[ind,2]];
        mm0 = Table[0,{sobolseqLength},{sobolseqLength}];
        mm1 = 3 IdentityMatrix[sobolseqLength];
        mm0[[;;polynomialDegree,;;polynomialDegree]] = 2;
        If[ highlight > 0,
            mm0[[;;, highlight;;highlight]] = 2;
            mm0[[highlight;;highlight, ;;]] = 2;
        ];
        ticks = {{{1,2,3,4,5,6,7,8,9,10,15,20,30},{1,2,3,4,5,6,7,8,9,10,15,20,30}}, {{1,2,3,4,5,6,7,8,9,10,15,20,30},{1,2,3,4,5,6,7,8,9,10,15,20,30}}};
        lbl = "sind"<>n2PaddedString[ind]<>StringJoin@(("_" <> ToString[#]) & /@ thisSetDirectionVectors);
        mmm = sobolMX + mm0 + mm1 ;
        MatrixPlot[mmm,FrameTicks->ticks,  ImageSize ->{sz,sz},
                ColorRules -> {0 -> White, 1 -> Blue, 2 -> (Lighter@Yellow), 3 -> Red,4 -> Green, 5 -> Green, 6 -> Green, 7 -> Magenta}, PlotLabel -> lbl, Mesh -> All ] 
    ]	(* getsobolMxGl*)

buildMSobolReturnMatrixPlot[indtab_,dbg_:True,modifiedDirectionVectors_:{},highlight_:0] :=
    Module[ {sobolDataFname,res},
        sobolDataFname = "data/sobol_init_tab.dat";  (* a copy of "data/new-joe-kuo-6.21201.dat";*)
        If[!FileExistsQ[sobolDataFname],Print[sobolDataFname," does not exist"];Abort[] ];
        sobolData = Drop[#,1]& @ Import[sobolDataFname];
        If[ modifiedDirectionVectors =!= {} && Length[indtab] != Length[modifiedDirectionVectors],
            Print["Length[indtab] != Length[modifiedDirectionVectors]"]; Abort[]
        ];
       res = Table[
	            If[ modifiedDirectionVectors === {},
	                getMsobol1DReturnMatrixPlot[ indtab[[i]], dbg,{},highlight],(*ELSE*)
	                getMsobol1DReturnMatrixPlot[ indtab[[i]],dbg,modifiedDirectionVectors[[i]],highlight ]
	            ]
        ,{i,Length[indtab]}];
        res
    ]

getMsobol1DReturnMatrixPlot[ind_:0,dbg_:False,modifiedDirectionVectors_:{},highlight_:-1] :=
    Module[ {j,msobol,thisSetDirectionVectors,sobolMX,polynomialDegree,polynomial},
        msobol = Table[1,{sobolseqLength}];
        polynomialDegree = sobolData[[ind,2]];
        polynomial = Join[{1}, IntegerDigits[#,2,polynomialDegree-1]& @ sobolData[[ind,3]], {1} ];

        If[modifiedDirectionVectors === {},
			thisSetDirectionVectors = sobolData[[ind,4;;]]; (* initialization of direction vectors as in new-joe-kuo-6.21201.txt http://web.maths.unsw.edu.au/~fkuo/sobol/ *)
        ,(*ELSE*)
			thisSetDirectionVectors = modifiedDirectionVectors;
        ];
        msobol[[;;polynomialDegree]] = thisSetDirectionVectors;
        Do[msobol[[i]] = BitOr[msobol[[i]],1],{i,sobolseqLength}];
        Do[
            msobol[[i]] = BitXor @@ Join[Table[2^(j) polynomial[[j+1]] msobol[[i-j]],{j,1,polynomialDegree}],{msobol[[i-polynomialDegree]]}];
        ,{i,polynomialDegree+1, sobolseqLength }];
        sobolMX = Table[0,{sobolseqLength},{sobolseqLength}];
        Do[ sobolMX[[i,;;i]] = (IntegerDigits[#,2,i]& @ msobol[[i]])  ,{i, sobolseqLength }];
        {sobolMX,msobol[[;;polynomialDegree]],getsobolMxGL[msobol,sobolMX,ind,thisSetDirectionVectors,highlight,300]}
    ] (* getMsobol1D *)


getMsobol1DdegreeN[degree_,ind_,dbg_:False,modifiedDirectionVectors_:{},highlight_:-1] :=
    Module[ {j,msobol,sobolM,sobolMX,mmm,mm0,mm1,mi2s[nSPP,7],thisSetDirectionVectors,ticks,lbl,polynomialDegree,data,fname},
        msobol = Table[1,{sobolseqLength}];
        mm1 = IdentityMatrix[sobolseqLength];
        polynomialDegree = degree;
        fname = "data/primitivePolynomialsDegree"<>ToString[degree]<>"_digits_sel100.dat";
        If[ !FileExistsQ[fname], Print[fname," does not exist"]; Abort[] ];
        data = Import[fname];
        polynomial = data[[Mod[ind,Length[data],1] ]];

        If[modifiedDirectionVectors === {},
			thisSetDirectionVectors = getRandomDirectionVectorsDegreeN[degree]; (* initialization of direction vectors as in new-joe-kuo-6.21201.txt http://web.maths.unsw.edu.au/~fkuo/sobol/ *)
        ,(*ELSE*)
			thisSetDirectionVectors = modifiedDirectionVectors;
        ];
        msobol[[;;polynomialDegree]] = thisSetDirectionVectors;
        (*If[ dbg,Print["getMsobol1D:  ind=",ind," degrees of primitive polynomials=",polynomial," of max degree=", polynomialDegree];
                Print["direction vectors:", msobol[[;;polynomialDegree]] ] ];*)
                mm2 = Table[0,{sobolseqLength},{sobolseqLength}];
                Do[ mm2[[i,;;i]] = (IntegerDigits[#,2,i]& @ msobol[[i]])  ,{i, sobolseqLength }];
        Do[msobol[[i]] = BitOr[msobol[[i]],1],{i,sobolseqLength}];
        Do[
            msobol[[i]] = BitXor @@ Join[Table[2^(j) polynomial[[j+1]] msobol[[i-j]],{j,1,polynomialDegree}],{msobol[[i-polynomialDegree]]}];
        ,{i,polynomialDegree+1, sobolseqLength }];
                sobolM = Table[0,{sobolseqLength},{sobolseqLength}];
                Do[ sobolM[[i,;;i]] = (IntegerDigits[#,2,i]& @ msobol[[i]])  ,{i, sobolseqLength }];
                (*sobolMX = T[sobolM];*)
                sobolMX = sobolM;
            If[dbg,
                        (*visumx1 = Table[ StringJoin[ToString /@ (Join[#,Table[" ",{sobolseqLength-i}]]& @ (IntegerDigits[#,2,i]& @ msobol[[i]]))],{i, sobolseqLength }];
                        visumx2 = Table[ StringJoin[ToString /@ (Join[Table[" ",{sobolseqLength-i}],#]& @ (IntegerDigits[#,2,i]& @ msobol[[i]]))],{i, sobolseqLength }];
                        Print[{Range[sobolseqLength]}//T//mf, visumx2//mf, visumx1//mf];
                        Print[mf @ sobolM];*)
                                mm0 = Table[0,{sobolseqLength},{sobolseqLength}];
                                mm0[[;;polynomialDegree,;;polynomialDegree]] = 2;
                                If[highlight > 0,
                                        mm0[[;;, highlight;;highlight]] = 2;
                                        mm0[[highlight;;highlight, ;;]] = 2;
                                ];
                        ticks = {{{1,2,3,4,5,6,7,8,9,10,16,24,32,40,48,56,64},None}, {None,{1,2,3,4,5,6,7,8,9,10,16,24,32,40,48,56,64}}};
                        lbl = "sind"<>n2PaddedString[ind-1]<>StringJoin@(("_" <> ToString[#]) & /@ thisSetDirectionVectors);
                        mmm = sobolM+mm0+mm1+mm2;
                        Print[MatrixPlot[mmm,FrameTicks->ticks,  ImageSize ->{ 600,600},
                                ColorRules -> {0 -> White, 1 -> Blue, 2 -> (Lighter@Yellow), 3 -> Green,4 -> Red, 5 -> Green, 6 -> Green}, PlotLabel -> lbl, Mesh -> All ] ];

        ];
        Return[{polynomial,msobol,sobolMX}]
    ] (* getMsobol1DdegreeN *)


vanderCorput[n_,base_:2,ndigits_:32] := With[{binaryCode=IntegerDigits[n,base,ndigits]}, FromDigits[Reverse@binaryCode,base]/base^ndigits ]
base2VDCInteger[n_,ndigits_:16]    :=   (FromDigits[#,2]& @ (Reverse @ (IntegerDigits[#,2,ndigits]& @ n)))

sobol1D[dim_,n_] := (* original Sobol's construction *)
    Module[ {i,seq = IntegerDigits[n,2] //Reverse,seqlen,nn},
        seqlen = Length[seq];
        If[dim == 0, (* Simple Inv-Radical *)
        	nn = BitXor @@ Table[ seq[[i]] 2^(seqlen-i), {i, seqlen}];
        	Return[ FromDigits[IntegerDigits[nn,2,seqlen] ,2] / 2^seqlen ]
        ,(*ELSE*) (* std Sobol *)
	        nn = BitXor @@ Table[ seq[[i]] msobol[[dim,i]] 2^(seqlen-i), {i, seqlen}];
	        Return[ FromDigits[IntegerDigits[nn,2,seqlen] ,2] / 2^seqlen ]
        ]
    ]

sobol1DInv[dim_,n_] := (* sobol1DInv[sobol1D[n]] == n *)
    Module[ {i,seq = IntegerDigits[n,2] ,seqlen,nn},
        seqlen = Length[seq];
        nn = BitXor @@ Table[ seq[[i]] msobol[[dim,i]] 2^(i-1), {i, seqlen}];
        Return[ FromDigits[Reverse@IntegerDigits[nn,2,seqlen] ,2] / 2^seqlen ]
    ]

sobol2d[n_] := {sobol1D[0,n],sobol1D[1,n]}

sobol1DGrayCode[n_] := sobol1DGrayCode[1,n]
sobol1DOriginal[n_] := sobol1D[1,n]

sobol2dGrayCode[n_] := sobol1DGrayCode[#,n]& /@ Range[0, 1]
sobol2dOriginal[n_] := sobol1D[#,n]& /@ Range[0, 1]

sobol3DGrayCode[n_] := sobol1DGrayCode[#,n]& /@ Range[0, 2]
sobol3DOriginal[n_] := sobol1D[#,n]& /@ Range[0, 2]

sobol4DGrayCode[n_] := sobol1DGrayCode[#,n]& /@ Range[0, 3]
sobol4DOriginal[n_] := sobol1D[#,n]& /@ Range[0, 3]


sobol5dGrayCode[n_] := sobol1DGrayCode[#,n]& /@ Range[0, 4]
sobol5dOriginal[n_] := sobol1D[#,n]& /@ Range[0, 4]


sobol6DGrayCode[n_] := sobol1DGrayCode[#,n]& /@ Range[0, 5]
sobol6DOriginal[n_] := sobol1D[#,n]& /@ Range[0, 5]

digits2str[digits_] := StringJoin[ToString /@ digits] 
str2n[str_] := FromDigits[#, 2] & @ Table[StringTake[str, {i}] // ToExpression, {i, StringLength[str]}]

getSobol1DPt[n_] := (* original Sobol's construction *)
    Module[ {i,seq = IntegerDigits[n,2] //Reverse,seqlen,nn},
        seqlen = Length[seq];
        nn = BitXor @@ Table[ seq[[i]] msobol[[1,i]] 2^(seqlen-i), {i, seqlen}];
        Return[ FromDigits[IntegerDigits[nn,2,seqlen] ,2] / 2^seqlen ]
    ]
    
getSobolPts2d[npts_,sobolIndTab_:{0,1}] := (* original Sobol's construction *)
    Module[ {},
    		{msobol,sobolMX} = buildMSobol[sobolIndTab,False];
        Table[ sobol2dOriginal[i] ,{i,0,npts-1}]//N
    ]
    
sobolScrambledPermuts2pts1D[permutTree_,sobolInd_:1,dbg_:False] :=
    Module[ {xresdigits,xdigitsReversed,ind,resBit,newpts,nbits,npts,pts},
    	{msobol,sobolMX} = buildMSobol[{sobolInd}];
        If[ !known[msobol], {msobol,sobolMX} = buildMSobol[{sobolInd,sobolInd}] ];
        nbits = Length[permutTree];
        npts = 2^nbits;
        pts = Round[2^(nbits) ( Table[ getSobol1DPt[i] //N,{i,0,npts-1}] ) ];
        xresdigits = Table[
            xdigitsReversed =  IntegerDigits[#,2,nbits]& @ pts[[ii]];
            Table[
                ind = 1 + FromDigits[#,2]& @ xdigitsReversed[[;;ilevel-1]];
                resBit = BitXor @@ {permutTree[[ilevel,ind]], xdigitsReversed[[ilevel]]};
                resBit
            ,{ilevel,Length[xdigitsReversed]}]
        ,{ii,npts}];
        newpts = FromDigits[#,2]& /@ xresdigits;
        If[ dbg,
            {permutTree -> newpts}//Print;
        ];
        newpts
    ] (* sobolScrambledPermuts2pts1D *)
   
applyOwenPermutTree1DReversed[permutTree_,pts1D_,dbg_:True] :=
    Module[ {xresdigits,xdigitsReversed,ind,resBit,newpts,nbits,npts},
        nbits = Length[permutTree];
        npts = 2^nbits;
        If[npts != Length[pts1D], Print["applyOwenPermutTree1D: permutTree and pts1D are of enequal size:",{npts, Length[pts1D]}]; Abort[] ];
        xresdigits = Table[
            xdigitsReversed =  Reverse @ (IntegerDigits[#,2,nbits]& @ pts1D[[ii]]);
            Table[
                ind = 1 + FromDigits[#,2]& @ xdigitsReversed[[;;ilevel-1]];
                resBit = BitXor @@ {permutTree[[ilevel,ind]], xdigitsReversed[[ilevel]]};
                resBit
            ,{ilevel,Length[xdigitsReversed]}]
        ,{ii,npts}];
        newpts = FromDigits[#,2]& /@ (Reverse @ xresdigits);
        If[ dbg, {permutTree -> newpts}//Print ];
        newpts
    ] (* applyOwenPermutTree1D *)

sobolScrambledInteger2pts1D[npts_,sobolInd_:1]:= sobolScrambledPermuts2pts1D[getRandomPermutTree[Log[2,npts]],sobolInd] 
sobolScrambledInteger2pts2d[npts_,{sobolInd1_,sobolInd2_}]:= {sobolScrambledInteger2pts1D[npts,sobolInd1],sobolScrambledInteger2pts1D[npts,sobolInd2]}//T


(*------------------------- end of supprot for stdSobol -------------------------*)


