TITLE 'Stage 1 - 1DTD'
COORDINATES
	Cartesian1
VARIABLES
p (threshold=0.001),l (threshold=0.001),m (threshold=0.001),q (threshold=0.001),h (threshold=0.001),N (threshold=0.001)
SELECT
ERRLIM=1e-3
threads=1
DEFINITIONS
	Dp=1.0e6					{Dp}
	Mup=1.0e6				{Mup}
	d_p=1.0e3					{dp}

	Dq=1.0e6					{Dq}
	Muq=1.0e6				{Muq}
	d_q=1.0e3					{dq}
	
	Dl=1.0e4					{Dl}
	Mul=1.0e5					{Mul}
	d_l=1.0e2					{dl}

	Dm=1.0e2					{Dm}
	Chim=1.0e4				{Chim}	
	Mum=1.0e2				{Mum}
	NuN=1.0e1{------------------}{1.0e1}
	{NuN2=1.0e1}	
	d_m=1.0e0				{dm}

	Dh=1.0e4					
	Nuh=1.0e4{------------------}{1.0e4}			
	d_h=1.0e2
	Kappa=1.0e0

	DN=1.0e-2
	
	Sigmap1=1.0e5		{Sigmap}
	Sigmap2=1.0e4		{Sigmap}
	Betap=1.0e0				{Betap}

	Sigmaq=1.0e0

	Sigmal=1.0e3{if t<6.7 then 3.0e3 else 1.0e3}			{Sigmal}
	Alphal=1.0e1{---------------}{1.0e1}
	Gammal=1.0e-1	

	Sigmam=4.0e-3		{Sigmam}
	P0=1.0e-1					{P0}
 	A=1.0e-1					{A0}
	Alpham=1.0e1{---------------}{1.0e1}
	Gammam=1.0e-1

	Sigmah=8.0e2{if t<6.7 then 3.2e2 else 5.5e2}

	{LEN=100}

EQUATIONS
	p:		Dp*del2(p)=-Mup*l*m/(1.0+l)+d_p*p+dt(p)

	q:		Dq*del2(q)=-Muq*l*m/(1.0+l)+d_q*q+dt(q)
	
	l:		Dl*del2(l)=Mul*l*m/(1.0+l)+d_l*l+dt(l)

	h:		Dh*del2(h)=Nuh*h*N/(Kappa+h)+d_h*h+dt(h)
	
	m:	Dm*del2(m)=Chim*div(m*grad(l))+Mum*l*m/(1.0+l)-NuN*h*N/(Kappa+h)+d_m*m+dt(m)

	N:		DN*del2(N)=-Mum*l*m/(1.0+l)+NuN*h*N/(Kappa+h)+dt(N)

BOUNDARIES
	REGION 1
		{MESH_DENSITY=0.1}
		START(0)
			POINT LOAD(p)=Sigmap1*l/(Betap+l)+Sigmap2*q
			POINT LOAD(q)=-Sigmaq*q
			POINT LOAD(l)=Sigmal*(1.0+h/Alphal)/(1.0+h/Gammal)
			POINT LOAD(h)=Sigmah
			POINT LOAD(m)=IF (p<P0) THEN 0 ELSE Sigmam*(1.0+h/Alpham)/(1.0+h/Gammam)*(1.0+A*q)*(p-P0)
			POINT LOAD(N)=0
		LINE TO(1)
			POINT LOAD(p)=0
			POINT LOAD(q)=0
			POINT LOAD(l)=0
			POINT LOAD(h)=0
			POINT LOAD(m)=0
			POINT LOAD(N)=0

TIME 0 TO 15
PLOTS
	FOR T=0 BY 1 TO 15
		ELEVATION(p) FROM (0) TO (1)
		{ELEVATION(dx(p)) FROM (0) TO (1)}
		ELEVATION(q) FROM (0) TO (1)
		{ELEVATION(dx(q)) FROM (0) TO (1)}
		ELEVATION(l) FROM (0) TO (1)
		{ELEVATION(dx(l)) FROM (0) TO (1)}
		ELEVATION(m) FROM (0) TO (1)
		{ELEVATION(dx(m)) FROM (0) TO (1)}
		ELEVATION(h) FROM (0) TO (1)
		{ELEVATION(dx(h)) FROM (0) TO (0.1)}
		ELEVATION(N) FROM (0) TO (1)
		{ELEVATION(dx(N)) FROM (0) TO (0.1)}

		ELEVATION(p,dx(p),q,dx(q),l,dx(l),m,dx(m),h,dx(h),N,dx(N)) from (0) to (1) export format "#x,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12" file "data_2"

histories
history(p) at (0) {export format "#t#r,#i" file "l(0)"}
history(integral(p)) {export format "#t#r,#i" file "int(p)3"}
history(q) at (0) {export format "#t#r,#i" file "m(0)"}
history(integral(q)) {export format "#t#r,#i" file "int(q)3"}
history(l) at (0) {export format "#t#r,#i" file "p(0)"}
history(integral(l)) {export format "#t#r,#i" file "int(l)3"}
history(h) at (0) {export format "#t#r,#i" file "h(0)2"}
history(integral(h)) {export format "#t#r,#i" file "int(h)3"}
history(m) at (0) {export format "#t#r,#i" file "m(0)"}
history(integral(m)) {export format "#t#r,#i" file "int(m)"}
history(N) at (0) {export format "#t#r,#i" file "p(0)"}
history(integral(N)) {export format "#t#r,#i" file "int(N)_Fieg_Sigma_l"}
END