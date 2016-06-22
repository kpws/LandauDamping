(* ::Package:: *)

tmp=Get["/data1/petter/angularCoefficientsEuclidean_ks=4000_acc=1000_nphi=200_phimin=-0.0490874_phimax=0.549779"];
acc=1000;
\[Phi]s=tmp[[1]];
f=OpenWrite["/data1/petter/smallPhi",PageWidth->Infinity];
Write[f,acc];
Write[f,Length@tmp[[1]]];
Write[f,Length@tmp[[2,1]]];
Write[f,N@#]&/@\[Phi]s;
WriteLine[f,#]&/@ToString/@CForm/@Flatten[ReIm[Flatten[tmp[[2]]]]];
Close[f];
