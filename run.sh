#!/bin/bash 
#PBS -N 150Nd
#PBS -o run150Nd.out
#PBS -e run150Nd.err 
#PBS -l select=1:ncpus=1:host=cn1   
ELE=Al # ZKr #Ca #Mo # Ti # Cr # Fe # Ni # Zn #Ge # Se # Ti62 # Cr64 # Fe66 # Ni68 # Zn70 # # Ge72 # Se74 # Kr #Sr # Zr # Mg # Si # Si # Mg
A=36
Nf=08

BlockLeveln=00
BlockLevelp=00
bType=2
qusiLevel=${BlockLeveln}
NUC=${A}${ELE}


pathwork=~/dio_f90
cd $pathwork
outputCopyPath=${pathwork}/data
pathexe=${pathwork}/exe
mkdir $pathexe
cp $pathwork/exe/dio $pathexe/dio
mkdir -p ${pathexe}/input
mkdir -p ${pathexe}/output
rm -rf ${pathexe}/input/*
rm -rf ${pathexe}/output/*
#--------------------------------------------------------------------
#      generate input files: dic.dat, gcm.dat, betgam.dat
#--------------------------------------------------------------------
# parameter sets that could be used in the beyond mean-field calculations
# Force    =  PC,PC-F1               ! Parameterset of the Lagrangian
#  V0      = 
# Force    =  PC,PC-PK1              ! Parameterset of the Lagrangian
#  V0      =  349.500   330.000      ! pairing strength for delta-force
#  V0      =  314.550   346.500      ! 
#--------------------------------------------------------------------
cd ${pathexe}/input
echo create input file ....
cat <<EOF > dio.dat 
n0f,n0b  = ${Nf}  8                    ! number of oscillator shells
b0       = -2.448                   ! oscillator parameter (fm**-1) of basis
beta0    =  0.00                    ! deformation parameter of basis
betas    =  0.50                    ! deformation beta2 of W.S. potential
bet3s    =  0.00                    ! deformation beta3 of W.S. potential
maxi     =  400                     ! maximal number of iterations
xmix     =  0.50                    ! mixing parameter
inin     =  1                       ! 1 (calc. from beginning); 0 (read saved pot.) 
${ELE} $A                           ! nucleus under consideration
Ide      =  4  4                    ! Pairing control: 1. no  2. Frozen  3. G   4. Delta
Delta    =  0.000000  0.000000      ! Frozen Gaps (neutrons and protons)
Ga       =  0.000000  0.000000      ! Pairing-Constants GG = GA/A
Delta-0  =  2.000000  2.000000      ! Initial values for the Gaps
Vpair    =  349.500   330.000       ! pairing strength for delta force
Force    =  PC-PK1
icstr    =  2                       ! Quadratic constraint (no 0; beta2 1; b2+b3 2)
cspr     =  10.00                   ! Spring constant
cmax     =  1.000                   ! cutoff for dE/db
iRHB     =  0                       ! 0: BCS; 1: RHB
iBlock   =  2                       ! 0: Non-blocking; 1: Block the given energy level; 2: Block according to K^\pi
bln      =  ${BlockLeveln}          ! block level of Neutron (valid when iBlock=1)
blp      =  ${BlockLevelp}          ! block level of Proton  (valid when iBlock=1)
Kn       =  1                       ! K of block neutron, 1: 1/2, 2: 3/2 ...  (valid when iBlock=2)
Pin      =  -1                       ! Pi of block neutron, 1: + , -1: -       (valid when iBlock=2)
Kp       =  2                       ! K of block proton,  1: 1/2, 2: 3/2 ...  (valid when iBlock=2)
Pip      =  1                       ! Pi of block proton, 1: + , -1: -        (valid when iBlock=2)
bMethod  =  ${bType}                ! block method, 1: blocking -> convergence;   2:convergence -> block; 3: convergence -> block -> convergence
qsn      =  ${qusiLevel}            ! one quasiparticle state (Neutron) of odd-mass nuclei
crankCase=  2                       ! 1:Belyaev formula; 2: Nilsson formula; 3: Odd A formula
c-------------------------------------------------------------------
EOF
#-----------------------------------
# mesh points in deformation q-space
####################################  
cat <<EOF > b23.dat
    0.30  0.00   0.00
EOF


cd ../
echo -e "\033[32m run dio...\033[0m"
./dio

touch ./output/HFB.wfs

echo calculation is finished !

dataDir=${outputCopyPath}/${NUC}_B${BlockLeveln}_${BlockLevelp}_${Nf}
mkdir $dataDir -p
cp ./output/* $dataDir -rf
echo Done!
