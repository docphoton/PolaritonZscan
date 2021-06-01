##############
#   zscan4level_all_lens.m
#
#  four-level model for non-linear optical effects for inclusion
#      in a structured optic in thin-film zscan. Open and Closed zscan
# 
# vector is vv=(p0'0', p-1-1, p00, p11, p0'-1, p0'0, p0'1, p-10, p-11, p01, * ) Here 0, 0' and -1 are exctied states. 1 is teh ground state. We mostly ignore 0', maybe casting it far away in detuning''. 
#
#  For RSA we think about the 2 photon excitation as leading through 1 (ground state) to '0' to -1. so we are thinking naively as the -1 state as an opical interval above the 0' state, no direct optical transition between the -1 and 1, but relaxation between -1 down to 1 and not back! Also from -1 to 0 and not back. So have to implement decay appropriately. 
#      
#  Code now gets SA and RSA in different regimes of couplings. The couplings are controlled through the Clebsch's' c10, c10p (for ground to 1st excited state manifold) and   cm10, cm10p from the excited' state manifold 0, 0' to the top state -1. Turn the first pair on to get SA and, subsequently, the second pair on to get RSA. 
# Another way: simply change G from smaller than 01 decoherence (RSA) to larger than 01 decoherence (SA)
#
# runstring: octave zscan4level_all.m 1e9 550 2.7 98 file1 30 
#               ""            light interval   clebsc dyelayThick fileout angle
############## CHANGELOG: 
#      11/17/2020: Included rest of the optical geometry of the zcsan setup, model closed zscan as non-linear lensing event
#      2/20/2019 : changing it to compute closed Zscan as well
#      1/07/2019 : adapting for zscan of polariton assembly
#      5/10/2016 : Mods to make estimate of central intensity more realistic while absorbing. Also updated with DBR mirrors having somewhat random layer thicknesses. Finally added some simple naive smoothing at the end for more realistic graphs...inhomog broadening of signal.  
##############
clear
# SECTION I: input data:
n_uno = str2num(argv(){1});    # Input intensity max at z=0, or scan parameter
p_due = str2num(argv(){2});    # zscan illumination wavelength, in nm
r_tre = str2num(argv(){3});    # mixing rate in 0-0' excited state manifold', in GHz 
q_quatro =  str2num(argv(){4});# Dye layer thickness , in nm
outfilename = argv(){5};       # specify output file name 
c_cinque = str2num(argv(){6}); # angle in degrees
    ##### interrogation specs
pi = 4.0*atan(1.0);            # indiana's 22/7th'
phase_offset = 0.0;  
I0     = n_uno;                 # Average input laser power (W)
wstart = (n_uno-50)*250;       # +-1800 loop range, could be wavelength, could be distance , could be angle
wend   = (n_uno-50)*250;
wstart = 00.0; 
wend = wstart;     
wav    = p_due;                # wavelength in nm; 
num    = 1;                   # number of points in [wstart:wend] grid
conv   = 3.0*413.5;            # convert 1/nm into eV
angle  = c_cinque*pi/180.0;    # angle between sample and optical axis
mode   = 1;                    # polarization=1 for Spol (TE), -1 for Ppol (TM)
closedRings = 400;             # number of r^2 rings to include in the closed zscan: have to keep this fairly high else get curvature of the Iout as a function of z...so yeah, check this at end by doing multiple runs at different number of integration rings
brightRing = 15; 	             # ring to single out and make realy bright...fro aperture defect testing. 
toll   = 1.0e-2;               # tollerance to converge the internal annealing of intensity to.
toll1  = 1.0e-5;               # we take things this small and larger to be consistent with 0...out numerical stability of the null() command 
anneal = 0.80;                 # annealing parameter for convergence of the intensity in the center of the stucture.
    ##### zscan parameters
eccentricity = 0.0;            # distance mismatch between the center of the rotary stage and the optical axis
zr = 2000.0;                   # Rayleigh range of optical system, in microns = \pi w_0^2/lambda, spot size is naively w = w_0 sqrt(1+z^2/z_r^2), here based on a spot radius of 20 microns
w0 = (zr*wav/1000.0/pi)**0.5;  # waist radius equation...in microns
laserPowerConv =  1.0/1000.0/200e-15;  # to convert measured input power (W) peak watts. Uses fact that the Clark laser is running at 1KHz, and has a 200 fS pulsewidth
    ##### optical cavity 
xc = q_quatro;                 # thickness (nm) of dye layers
xm1 = 18.0;                    # thickness (nm) of front mirror
xm2 = 24.0;                    #    " " back mirror
n1_0=1.58+.000002*i;           # PS index of refraction
#n1_0=1.8*700.0/660.0;  	# TEMPORARY!!! To put a 134nm stack in same resonance position as the polariton example
lin_dispn1 =-6.71e-5;          #   linear dispersion of PS away from 566nm
n2_0=1.49+0.00*i;              # index of PMMA
lin_dispn2 =-3.07e-5;          #   linear dispersion of PMMA away from 554nm
  n1 = n1_0+lin_dispn1*(wav-566.0);    # naive simple approximation of dispersion over range of interest...visible spectrum. here and next line for plastics we're not using' 
  n2 = n2_0+lin_dispn2*(wav-554.0);
n3=1.58-.08i;                  # index of TTY 
nm_r = -0.01;                  # index of mirror material (Ag) from  McPeak, et. al.,  2015, 0.3-1.7 microns,
nm_r_slope = 0.0001;           # https://refractiveindex.info/?shelf=main&book=Ag&page=McPeak   
nm_i = -0.81;
nm_i_slope = +.008;            # + sign overall is absorption, - sign overall is gain
ns = 1.5;                      # glass index: the substrate the polariton is built on 
    #####  Center (the dye) layer parameters. 4 level model  
dens = 1.8e21*1e6;             # dye number density from Liu PhysRevB.92.155301 (2015), here per m^3
dipoleM = 7.3*3.336e-30/(1.054e-34*1e9);      # averaged debye matrix element of center layer. all from Liu PhysRevB.92.155301 (2015). Put in conversion factor (hbar*e) to convert MKS to GHz in Resulting Rabi frequency. This just totals out to 4320 GHz/(V/m)
nc0 = real(n1_0);              # base index (carrier material, off resonant response) ... dissolved in PMMA
dispersion = 2500;             # natural, linear dispersion in the base material, see below... 
wavec = 575; #620;             # center wavelength (nm) of absorption of the |1> -> |0> (exciton) from Liu PhysRevB.92.155301 (2015)
wavec2 = 3200;                 # center wavelength of the 0-> -1 absorption from the excited state
br0 = .3;                      # branching ratio of population decay from -1 into the 0 versus the 1 state (unique ground state) 	
br0p = .3;                     # branching ration pop decay from -1 into 0' state' 	
split0 = -94000.0;             # Splitting of the excited state manifold, in GHz, between the 0 and 0' states'			
delta = -2.0*pi*(1.0/wavec-1.0/wavec2)*3.0e8;       # two photon detuning between 1 -> -1 , in GHz, so note that  D0-delta is the detuning from the 0 to the -1 (i.e., between excited states). 
G  = 1.0e4;                    # -1 population decay rate, GHz						      
mix = r_tre;                  # how fast is mixing in the excited state? (GHz)
gam0 = 1.5e4;                  # pop decay rate for 0 level exctied state,
gam0p = 0.25e4;                # pop decay rate for 0-prime excited state
gam20 =  .5*(G+mix+gam0)+3.0e4;    # relaxation of coherence between 0 and -1 (excited states) 
gam2p0 = .5*(G+mix+gam0p)+7.0e4;   #  for the decoherence rate between 0 and 1 (groudn to 1st excited state) , last figure is the 0.29eV measured width in Liu PhysRevB.92.155301 (2015)
gam2e = mix+.5*(gam0+gam0p)+1.0e3; # coherences lifetime inverse between excited states
G2  = 0.5*G+1.50e6;            # -1 <-> 1 coherence width
cm10 = 0.0;                    # Clebsch-Gordon coefficients connecting up the reduced matrix element to the dipole for that state. 0 is the upper level at gap+split0, that is the 2F' = 2' level. We toggle this between 2.7 for RSA and .1 for SA! 
cm10p = 0.0*cm10;              # 0p is the lower level , just at gap energy, thatis the F' = 1' level , here a cluge for non-participating levels as compared to teh F'=2 final state' 
c10  =  2.0;                   # note could have sign differences, etc.  between the excited states here!  
c10p  = 0.01*c10; 
thickness = 0.2;               # optical thickness parameter...part of the temperature dependence in the data perhaps. 
n21 = 0.0; #8e-11/50.0;              # just test values n21 and n22
n22 = -2e-10/50.0;
n44= 0.0; 
k1  = 0.0; #2e-8;
#####  Setup arrays/vars
I0 = I0*laserPowerConv*2.0/(8.85e-12*3e8)*dipoleM*dipoleM;    # this should be (E*matrixLmnt)^2 after you divide by the area 
dens = dens*(dipoleM*1.054e-34*1e9)^2/(6.626e-34*3e8*8.85e-12)*3e8;    # this is the coefficient as used in the determination of the contribution the dye makes to the index of refraction
being = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0] ; # this implements p11+p00+p0'0'+p-1-1= 1.0
v = rand ("state");            # for use in the random layer thicknesses 
#This returns a column vector v of length 625. Later, you can restore the random number generator to the state v using the form
#rand ("state", v);
for w=1:num                    # initialize arrays
	xx(w) = wstart + (wend-wstart)*(w-1)/num; 
	R(w) = 0.0;                  # final answer arrays. 
	F(w) = 0.0;
  J(w) = 0.0;  
endfor
k0 = 2.0*pi/wav;                # note that this one is in nm^-1, since layer thicknesses are too
for w=1:num                              ##### MAIN LOOP
   w2 = w0*w0*(1.0+xx(w)*xx(w)/zr/zr); 
   area = wav*zr*w2/1000.0/w0/w0*1.0e-12;#  throat area in  meters^2 = pi*w2 
#  area = 1.0;                           #  for wavelength scan only!!
   I1 = I0/area;                         #  for the wavelength scan? old line leftover. using I = 2*(Power of pulse)/(pi w0^2) 
#  k0 = 2*pi/xx(w);                      #  for wavelength scan only!!
   k = k0;                                # the wavevector (vacuum)
#  wav = xx(w);                          #  for wavelength scan only!! 
  D0= (k-2.0*pi/wavec)*3.0e8;            # detuning in GHz (really GRad/s radians per second) 
  ll=1;
  ddrr = w2*9.0/closedRings;             # so start loop at rr=0 actually.
  rr = -ddrr;                            # and initialize the parameters for the run
  Efar = 0;
  popDiff  = 0 ; 
  PLpop = 0; 
  Itot = 0; 
  lightOut = 0; 
  loopThrough = 0; 
while(loopThrough<closedRings+1)         # loop for integration needed for closed zscan
  rr = rr+ddrr;          # add up each r^2 ring in the closed Zscan integral 
  loopThrough = loopThrough+1;  
    I1 = 2.0*I0/pi/w2*exp(-2.0*rr/w2)*1.0e12;   #   get intensity in the rings we're' adding up.
    Ireset = I1;                         # for cavity update...see below in update step. 
  stop = 0; 
  while(stop==0)                         # iteration to get non-linear light field center approx right. 
    ###### 4 LEVEL CODE: put in the solution of the 4 level quantum optics model.
#while(0)					   ##TEMPORARY to remove the 4-level system computation
gam = gam0;                              # excited state pop relaxation parameters
gamp = gam0p;                            #  ''
gam2 = gam20;
gam2p = gam2p0;  
D = D0+split0-delta/2.0; 
omega = real(sqrt(I1)); 
A = real(sqrt(I1));
MM=zeros(16); #  1      2     3     4     5    6     7     8     9    10    11      12    13     14    15   16
#          vv=(p0'0', p-1-1, p00, p11, p0'-1, p0'0, p0'1, p-10, p-11, p01, * p-10', p00', p10', p0-1, p1-1, p10)	       
MM(1, :) = [-gamp-mix, G*br0p, mix, 0, i*cm10p*conj(A), 0, i*c10p*omega,0, 0,0, -i*cm10p*A, 0, -i*c10p*omega, 0, 0,0];          #
MM(2, :) = [0, -G, 0, 0, -i*cm10p*conj(A), 0, 0, i*cm10*A, 0, 0, i*cm10p*A, 0, 0, -i*cm10*conj(A), 0, 0];       #
MM(3, :) = [mix, G*br0, -gam-mix, 0, 0, 0, 0,  -i*cm10*A, 0,  i*c10*omega, 0, 0, 0, +i*cm10*conj(A), 0, -i*c10*omega ];        #
MM(4, :) = [gamp, G*(1.0-br0p-br0), gam, 0.0, 0, 0, -i*c10p*omega, 0, 0, -i*c10*omega, 0, 0, i*c10p*omega, 0, 0, i*c10*omega];        #
MM(5, :) = [i*A*cm10p, -i*A*cm10p, 0, 0, i*D-gam2p-i*delta/2, i*cm10*A, 0, 0, 0, 0, 0, 0, 0, 0, -i*c10p*omega, 0];         #       
MM(6, :) = [0, 0, 0, 0, i*cm10*conj(A), i*split0-gam2e, i*c10*omega, -i*cm10p*A, 0, 0, 0, 0, 0, 0, 0, -i*c10p*omega];    #              
MM(7, :) = [i*c10p*omega, 0, 0, -i*c10p*omega, 0, i*c10*omega, i*D-gam2p+i*delta/2, 0, -i*cm10p*A, 0, 0, 0, 0, 0, 0, 0];   #         
MM(8, :) = [0, i*cm10*conj(A), -i*cm10*conj(A), 0, 0, -i*cm10p*conj(A), 0, -gam2+i*delta/2-i*(D-split0), i*c10*omega, 0, 0, 0, 0, 0, 0, 0];     # 
MM(9, :) = [0, 0, 0, 0, 0, 0, -i*cm10p*conj(A), i*c10*omega, i*delta-G2, -i*cm10*conj(A), i*c10p*omega, 0, 0, 0, 0, 0];     #       
MM(10, :)= [0, 0, i*c10*omega, -i*c10*omega, 0, 0, 0, 0, -i*cm10*A, i*(D-split0)+i*delta/2-gam2, 0, i*c10p*omega, 0, 0, 0, 0];       # 
MM(11, :)= [-i*cm10p*conj(A), i*cm10p*conj(A), 0, 0, 0, 0, 0, 0, i*c10p*omega, 0, -i*D-gam2p+i*delta/2, -i*cm10*conj(A), 0, 0, 0, 0];   #
MM(12, :)= [0, 0, 0, 0, 0, 0, 0, 0, 0, i*c10p*omega, -i*cm10*A, -i*split0-gam2e, -i*c10*omega, i*cm10p*conj(A), 0, 0];        #
MM(13, :)= [-i*c10p*omega, 0, 0, i*c10p*omega, 0, 0, 0, 0, 0, 0, 0, -i*c10*omega, -i*D-gam2p-i*delta/2, 0, i*cm10p*conj(A), 0];  #
MM(14, :)= [0, -i*cm10*A, i*cm10*A, 0, 0, 0, 0, 0, 0, 0, 0, i*cm10p*A, 0, i*(D-split0)-gam2-i*delta/2, -i*c10*omega, 0];   #
MM(15, :)= [0, 0, 0, 0, -i*c10p*omega, 0, 0, 0, 0, 0, 0, 0, i*cm10p*A, -i*c10*omega, -i*delta-G2, i*cm10*A];                #
MM(16, :)= [0, 0, -i*c10*omega, i*c10*omega, 0, -i*c10p*omega, 0, 0, 0, 0, 0, 0, 0, 0, i*cm10*conj(A), -i*(D-split0)-i*delta/2-gam2];   #

vv=null(MM,1);                            # this only makes sense since physically there can only be a single null vector.
scale = being*vv; 
vv = vv./scale;
    ###### END 4 LEVEL CODE
    ###### Sanity checks...are all populations positive? MASK lights up problems! 
if (abs(imag(vv(1))>toll1)||(abs(imag(vv(2)))>toll1)||(abs(imag(vv(3)))>toll1)||(abs(imag(vv(4)))>toll1))
 mask3=0.0
 else
   mask3=1.0; 
endif
if ((real(vv(1))<-toll1)||(real(vv(2))<-toll1)||(real(vv(3))<-toll1)||(real(vv(4))<-toll1))
  mask1=0.0                               # error....blanking mask
  vv(1)
  vv(2) 
  vv(3) 
  vv(4) 
 else
   mask1=1.0; 
vv(1) = abs(vv(1)); 
vv(2) = abs(vv(2));
vv(3) = abs(vv(3)); 
vv(4) = abs(vv(4)); 
endif 
for ll=5:10                               # check conjugation symmetry
	 if ((abs(real(vv(ll))-real(vv(ll+6)))+abs(imag(vv(ll))+imag(vv(ll+6))))>toll1)
	   mask4=ll
             real(vv(ll))
	           real(vv(ll+6))
             imag(vv(ll))
             imag(vv(ll+6))
	 else
	   mask4=1.0; 
         endif
endfor
if(abs(vv(5)*conj(vv(5)))-abs(vv(1)*vv(2))>toll1)  # check schwartz identity...six of them
  mask5=0.0
  else
    mask5=1.0; 
endif
if(abs(vv(6)*conj(vv(6)))-abs(vv(1)*vv(3))>toll1)
  mask6=0.0
  else
    mask6=1.0; 
endif
if(abs(vv(7)*conj(vv(7)))-abs(vv(1)*vv(4))>toll1)
  mask7=0.0
  else
    mask7=1.0; 
endif
if((abs(vv(8)*conj(vv(8)))-abs(vv(3)*vv(2)))>toll1) 
  mask8=0.0
  else
    mask8=1.0; 
endif
if((abs(vv(9)*conj(vv(9)))-abs(vv(4)*vv(2)))>toll1) 
  mask9=0.0
  else
    mask9=1.0; 
endif
if((abs(vv(10)*conj(vv(10)))-abs(vv(3)*vv(4)))>toll1) 
  mask10=0.0
  else
    mask10=1.0; 
endif
    ###### DONE with sanity checks.... 
# INDEX: Now compute the atomic/molecular polarization's contribution to the (nonlinear) index of refraction'  
#    Note the check of the units here. The contribution to the dielectric function is dimeless. The term is (n d^2/(eps_0\hbar)*\rho_{off}\omega ! Works! 
n = nc0*nc0-dens*1.0e-9*(c10*vv(10)/omega+c10p*vv(7)/omega+cm10*vv(8)/A+cm10p*vv(11)/A);   # square it first...will square root it in dielectric function.
# END INDEX 
if(real(n)<0.0)
   n=.0000000001;                                    # error if go to far with pumping up response , with dens or with I0! This flags region. 
   mask2=0.0;                                        # output blanking mask
 else
   mask2=1.0; 
endif
  nc = sqrt(n);
  if (loopThrough==1)
    nc_center = real(nc); 
  endif
#endwhile                                            #TEMPORARY: REmove 3-level polariton test
#nc = n1_0+i*(k1+n21*I1)+n22*I1+n44*I1*I1;            # for testing only.
#        mask1=1.0;
#        mask2 = 1.0; 
  cos_theta = sqrt(nc*nc-sin(angle)**2)/nc;          # 
  mult=nc*cos_theta**mode;            
  k = k0*cos_theta; 
  c1=cos(nc*k*xc);
  s1=sin(nc*k*xc);                                   #
  center=[c1,i*s1/mult;i*s1*mult,c1];    
nm = nm_r+nm_r_slope*wav+i*(nm_i+nm_i_slope*wav);    # index model for the Ag mirrors. 
nm1 = nm; 
cos_theta = sqrt(nm1*nm1-sin(angle)**2)/nm1;         # this method of doing transfer matrix for a medium with a complex refractive index at different incidence angles is due to Petr Kuzel...nice!! 
mult=nm1*cos_theta**mode;            
k = k0*cos_theta; 
c1=cos(nm1*k*xm1);
s1=sin(nm1*k*xm1);                                   # first mirror done now do second
mm1=[c1,i*s1/mult;i*s1*mult,c1]; 
nm2 = nm; 
cos_theta = sqrt(nm2*nm2-sin(angle)**2)/nm2;           
mult=nm2*cos_theta**mode;            
k = k0*cos_theta; 
c1=cos(nm2*k*xm2);
s1=sin(nm2*k*xm2);                                   # second mirror
mm2=[c1,i*s1/mult;i*s1*mult,c1]; 
  mtot = mm2*center*mm1 ;    # this is the full transfer matrix of the system
# done with full mtot: now compute observables, for mtot on substrate
 mult0 = 1.0*cos(angle)**mode;                      # air coefs
  cos_theta = sqrt(ns*ns-sin(angle)**2)/ns;          # substrate coefs
  mult_final=ns*cos_theta**mode;
  trans  = 2.0*mult0/(mm2(1)*mult_final+mm2(4)*mult0-mm2(3)*mult0*mult_final-mm2(2));  # last mirror transmissitivity
  TT  = 2.0/(mtot(1)*mult_final+mtot(4)*mult0-mtot(3)*mult0*mult_final-mtot(2));       # whole system transmissitivity
  RR  = (mtot(1)*mult_final-mtot(4)*mult0+mtot(3)*mult0*mult_final-mtot(2))*TT/2.0/mult0; 
  cos_theta = sqrt(nc*nc-sin(angle)**2)/nc;   
  k = k0*cos_theta;
   Iprime = Ireset*real(TT*conj(TT))/real(trans*conj(trans))*cos(angle)/real(cos_theta);   % estimate of the intensity just inside end of the cavity in that ring 
   Iprime = Iprime*(1.0+exp(2.0*imag(nc)*k*xc))/2.0;                                % average it out with the intensity in front of cavity in that ring
     ccatch = abs(Iprime-I1)/I1;
   if (ccatch<toll) % test to see if we've' converged on a self-consistent approximate solution
       stop=1;                                      # Done! Converged!
     else
       I1 = I1+anneal*(Iprime-I1); 
       stop=0;                                      # not done...keep iterating the intensity
   endif
   #stop = 1;                                       # temporarily, if do not want to itterate while testing 
 endwhile                                           # END while for internal field determination 
 popDiff = popDiff+(v(4)-v(3))*ddrr/v(4); 
 PLpop = PLpop+(vv(1)*vv(4)-abs(vv(7))*abs(vv(7)))*ddrr/vv(4)/vv(4); 
 #PLpop = PLpop+vv(1)*ddr/vv(4); 
aperture = 1.0; 
if(loopThrough==brightRing)
    aperture=1.0;
endif
 Itot = Itot+I1*ddrr*aperture; 
 lightOut = lightOut+I1*abs(TT)**2*ddrr*aperture; 
 Efar = Efar + aperture*ddrr*TT*exp( (-1.0/w2+i*pi*1000.0/wav*w0*w0*xx(w)/zr/zr/w2)*rr);
 lightThrough(loopThrough) = I1*abs(TT)**2*ddrr;
 location(loopThrough) = (rr+ddrr)**0.5; 
 phaser(loopThrough) = -atan(imag(TT)/real(TT))+phase_offset;
 if((loopThrough>2)&&(phaser(loopThrough)-phaser(loopThrough-1)>pi/3.0)) # have to unwind the phase of the naive atan
  phase_offset = phase_offset-pi;
  phaser(loopThrough) = phaser(loopThrough)-pi; 
 endif 
 if((loopThrough>2)&&(phaser(loopThrough)-phaser(loopThrough-1)<-pi/3.0)) # have to unwind the phase of the naive atan
  phase_offset = phase_offset+pi;
  phaser(loopThrough) = phaser(loopThrough)+pi; 
 endif 
endwhile                                            # END while for closed zscan r^2 integral. 
EE = I0*abs(Efar*conj(Efar))/w2/closedRings/closedRings;     # final far field closed zscan INTENSITY
J(w) = (1.0-abs(TT)**2-abs(RR)**2)*mask1*mask2;     # FINAL REFACTORING total absorbed percentage
R(w) = lightOut*mask1*mask2/closedRings;            # total light: open zscan
F(w) = EE*mask1*mask2;                              # far field intensity, stored in array F
J(w) = nc_center;                                       # END FINAL REFACTORING, excited state population,                                                  
J(w) = R(1); 
#printLight = R(w)
  endfor                                            # END main loop for open zscan. Comment this out for closed/open zscan! See below!! 
xx = xx/zr;                                         
disp(real(Itot))
arraz = [real(xx); real(R); real(F); real(J)];      # REPORTING
#arraz = [location; phaser];
#plot (xx, R);
csvwrite (outfilename, arraz');  
# post processing section for focal lengths' 
[xmax, imax] = max (phaser);
[xmin, imin] = min (phaser);
rmcutoff = phaser(closedRings)-abs(imax-imin)*(xmax-xmin)/20.0/(imax-imin); 
for ww=1:closedRings
  ttemp(ww) = (phaser(ww)-rmcutoff)*(phaser(ww)-rmcutoff);  # use a new unused array so no confusion determining length, and what parts are min and max
endfor
[ymin, rmi] = min(ttemp); 
rm = location(rmi); 
lightinLens = 0.0; 
lightoutsideLens = 0.0; 
for ww=1:rmi
   lightinLens = lightinLens+lightThrough(ww); 
endfor
for ww=rmi:closedRings
   lightoutsideLens = lightoutsideLens+lightThrough(ww); 
endfor
fprime = rm*rm*k0*1000.0/(phaser(rmi)-phaser(1));     # the focal length in microns
#disp(fprime)
fprime = fprime/1000.0;         # microns to mm 
z = xx(w)/1000.0;               # microns to mm
f2 = 50.0;                      # focal length of collecting lens (mm)
mm = 120.0;                     # distance to closed channel aperture (mm)
closedOverOpen = (lightinLens/(1.0+(z*z*(1.0/f2-mm/f2/f2)-z)/fprime)**2+lightoutsideLens)/(lightinLens+lightoutsideLens)-1.0; 
#disp(fprime)
#disp(real(closedOverOpen))
#endfor                   #  this is end of loop for closed/open zscan, uncomment it for that. 
#disp(xx(imax)-xx(imin)) 
#csvwrite ("output", arraz');
#csvwrite ("state", v'); 
#axis([600,620,0,2]);
# plot(xx,Ra2,xx,Ta2,xx,Rb2,'p',xx,Tb2,'d');                  # END REPORTING

#  smooth = 1;           # for smoothing output stream if in a hurry and not averaging enough over ensemble
# naive smoothing to remove burrs from randomness...recognition of the inhomogeneous broadening effect of the transverse non-uniformity
#while(smooth != 1)
#  for w=smooth+1:num-smooth
#          AAA = 0;
#          RRR = 0; 
#          FFF = 0;
#    for jj=-smooth:smooth
#       AAA = AAA+A(w+jj);
#       RRR = RRR+R(w+jj);
#                   FFF = FFF+F(w+jj); 
#    endfor
#          B(w) = AAA/(2*smooth+1); 
#          S(w) = RRR/(2*smooth+1); 
#          GG(w) = FFF/(2*smooth+1); 
#  endfor
#  smooth=1; 
#endwhile  

