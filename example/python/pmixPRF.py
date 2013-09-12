#! /usr/bin/env python
def usagemenu():
    print 'syntax: ./pmix.py 1:ON 2:Phi [3:printtype]'
    print '        ON    -> octane number: percent (by volume) of iso-Octane in a mixture of'
    print '                 iso-Octane and n-Heptane'
    print '        Phi   -> equivalence ratio (0:10) for a mixture of fuel-air'
    print '                 fuel: mixture of iso-Octane and n-Heptane with ON defined above'
    print '        printtype -> (optional) 1-output mole fractions for ignition code input'
    print '                                                                     '
    print '        Air is assumed to be (N2:78.08%,Ar:0.97%,O2:20.95%) by volume'
    print '                                                        '
    print '                  1    3n+1      0.7808      0.0097     '
    print '        C H    + --- * ---- (O + ------ N  + ------ Ar) '
    print '         n 2n+2  phi    2     2  0.2095  2   0.2095     '
    print '                                                        '    
# retrieve input parameters and test correctness
import sys

if ( len(sys.argv) < 3 ) :
    usagemenu()
    sys.exit()

ON  = float(sys.argv[1])
phi = float(sys.argv[2])

onrat = ON/100

outtype = 0
if ( len(sys.argv) == 4 ) :
    outtype = int(sys.argv[3])

if  ( phi <= 0 ) | ( phi > 10) :
    print 'Error : illegal value for equivalence ratio : ',phi
    usagemenu()
    sys.exit()

# define element mass (taken from ckinterp)
ne = 5
EW = [1.00797, 15.99940, 12.01115, 14.00670, 39.94800] # H O C N AR

# define species compositions and names
iFiO = [2*8+2, 0, 8, 0, 0]
iFnH = [2*7+2, 0, 7, 0, 0]
iO2  = [0    , 2, 0, 0, 0]
iH2O = [2    , 1, 0, 0, 0]
iCO2 = [0    , 2, 1, 0, 0]
iN2  = [0    , 0, 0, 2, 0]
iAR  = [0    , 0, 0, 0, 1]

ispec = {'iO':iFiO, 'nH':iFnH, 'O2' :iO2 , \
         'CO2':iCO2, 'H2O':iH2O, \
         'N2' :iN2 , 'AR' :iAR}

# define N2/O2 mole ratio in air
#n2o2air = 0.79/0.21
n2o2air = 0.7808/0.2095
aro2air = 0.0097/0.2095

# compute molecular weights
for spec in ispec:
    sumW = 0.0
    i    = 0
    while i < ne:
        sumW += ispec[spec][i]*EW[i]
        i += 1
    ispec[spec].append(sumW)
Wfio = ispec['iO' ][ne]
Wfnh = ispec['nH' ][ne]
Wo2  = ispec['O2' ][ne]
Wn2  = ispec['N2' ][ne]
Wco2 = ispec['CO2'][ne]
Wh2o = ispec['H2O'][ne]
War  = ispec['AR' ][ne]

small  = 1.e-15

# compute REACTANTS
alpha  = (11+3*onrat/2) / phi 
betaN2 = alpha * n2o2air
betaAR = alpha * aro2air

denom   = Wfio*onrat+Wfnh*(1-onrat)+alpha*Wo2+betaN2*Wn2+betaAR*War
Yr_io = onrat    *Wfio/denom
Yr_nh = (1-onrat)*Wfnh/denom
Yr_o2 = alpha    *Wo2 /denom
Yr_n2 = betaN2   *Wn2 /denom
Yr_ar = betaAR   *War /denom

denom   = 1.0+alpha+betaN2+betaAR
Xr_io = onrat    /denom
Xr_nh = (1-onrat)/denom
Xr_o2 = alpha    /denom
Xr_n2 = betaN2   /denom
Xr_ar = betaAR   /denom

print 'spec IC8H18  ',Xr_io
print 'spec NC7H16  ',Xr_nh
print 'spec O2      ',Xr_o2
print 'spec N2      ',Xr_n2
print 'spec AR      ',Xr_ar

if outtype == 1:
    print 'Y_IC8H18  ',Yr_io
    print 'Y_NC7H16  ',Yr_nh
    print 'Y_O2      ',Yr_o2
    print 'Y_N2      ',Yr_n2
    print 'Y_AR      ',Yr_ar
    print 'sum of mole fractions:',Xr_io+Xr_nh+Xr_o2+Xr_n2+Xr_ar
    print 'sum of mass fractions:',Yr_io+Yr_nh+Yr_o2+Yr_n2+Yr_ar
