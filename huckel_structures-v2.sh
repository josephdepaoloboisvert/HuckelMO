#only print real part of eigenvalues   Feb 28th 2018 Jingbai Li
#V2-UI revised, interfaced with GView Aug 25th 2017 Jingbai Li
#Huckel_structures_calculator Feb 22th 2016 Jingbai Li 

#!/bin/sh
#Variables

#home=~/Desktop/Huckel_Structures
home=./
cd $home


#Title section
echo "     ===================================================================="
echo "     ||                                                                ||"
echo "     ||                                                                ||"
echo "     ||           HUCKEL Molecular Orbital (HMO) Calculations          ||"
echo "     ||                                                                ||"
echo "     ||                 CHEM344  Physical Chemistry II                 ||"
echo "     ||                    Department of Chemistry                     ||"
echo "     ||                Illinois Institute of Technology                ||"
echo "     ||                                                                ||"
echo "     ||                        By Jingbai Li                           ||"
echo "     ||                     jli129@hawk.iit.edu                        ||"
echo "     ||                                                                ||"
echo "     ||                 Version 2.0  Aug 30th 2017                     ||"
echo "     ||                                                                ||"
echo "     ===================================================================="
echo "                                                                              "
echo "=============================================================================="
files=`ls *.mol2|wc -l|awk '{print $1}'`
if [ "$files" = "0" ]
then
echo "                            NO INPUT FILES FOUND   "
exit
else
echo "                             SELECT INPUT FILES    "
fi
echo "=============================================================================="


for ((f=1;f<=$files;f++))
do
f[$f]=`ls *.mol2|sed -n "$f"p`
echo " $f. ${f[$f]}"
done
echo "=============================================================================="
read -e in
if [ "$in" -gt "$files" ]
then
echo "=============================================================================="
echo "                            NO INPUT FILES FOUND   "
echo "=============================================================================="
exit
fi
input="$home/${f[$in]}"
echo "=============================================================================="
echo "  GaussianView Molecular File:                                                "
echo "  $input                                                                      "
echo "=============================================================================="
echo "                                                                              "
echo "                                                                              "
echo "=============================================================================="
echo "                       BUILD UP SECULAR DETERMINANT                           "
echo "=============================================================================="
#

#Initial Parameters
alpha=1.292              #electrostatic
beta=-3.26934291633      #covalent
ele_num=0                #electrons
orb=0                    #orbitals

#alpha beta types: 1:C=C-C 2:C-N-C 3:C=N-C



#diagonal-alpha, alpha’=alpha + h * beta, C=C:0 C=N:2 C-N:1.5
#connected-beta, beta’=k*beta, C=C:1 C=N:1 C-N:0.8

atomsec=`sed -n '/<TRIPOS>ATOM/,/<TRIPOS>BOND/p' "$input"|sed 1d|sed '$d'`
atom=`echo "$atomsec"|wc -l|awk '{print $1}'`

bondsec=`sed -n '/<TRIPOS>BOND/,$p' "$input"|sed 1d|awk '{if ($4==1) $4="S";else if ($4==2) $4="D";else $4="N"; print $1, $2, $3, $4;}'`
bond=`echo "$bondsec"|wc -l|awk '{print $1}'`


#detect atom type, matrix dimension, diagonal element
for ((a=1;a<=${atom};a++))
do
  atom_type=`echo "$atomsec"|sed -n "$a"p|awk '{print $(NF)}'`  
  if [ "$atom_type" = "C" ]
  then 
  valence=`echo "$bondsec"|awk '{print $2, $3}'|grep -w "$a"|wc -l|awk '{print $1}'` #valence bond
    if [ "$valence" -lt "4" ]
    then
    typ[$a]=1   # type 1 is conjugated carbon atom
    h=0.0000
    ele_num=`expr $ele_num + 1`
    orb=`expr $orb + 1`
    else
    typ[$a]=0
    fi
  elif [ "$atom_type" = "N" ]
  then
  pyridine=`echo "$bondsec"|awk '{print $2, $3, $4}'|grep -w "$a"|grep -w "D"`
    if [ "$pyridine" ]
    then 
    typ[$a]=3   # type 3 is conjugated pyridine N atom (C=N-C)
    h=2.0000
    ele_num=`expr $ele_num + 1`
    orb=`expr $orb + 1`
    else
    typ[$a]=2   # type 2 is conjugated pyrolle N atom (C-N-C)  
    h=1.5000
    ele_num=`expr $ele_num + 2`
    orb=`expr $orb + 1`
    fi
  else
  typ[$a]=0     # type 0 is out of conjugation
  h=`echo "scale=10;(-1*$alpha/$beta)"|bc`
  fi
  m=`expr \( $a - 1 \) \* $atom + $a`
  matrix[$m]=`echo "scale=10;($alpha+$h*$beta)"|bc`
done  




#detect bond type, none-zero matrix element

for ((b=1;b<=${bond};b++))
do
  b1=`echo "$bondsec"|sed -n "$b"p|awk '{print $2}'`    #atom1
  b2=`echo "$bondsec"|sed -n "$b"p|awk '{print $3}'`    #atom2
  if [ "${typ[$b1]}" = "0" ] || [ "${typ[$b2]}" = "0" ]
  then
  k=0.0000     # type 0 is out of conjugation 
  elif [ "${typ[$b1]}" = "2" ] || [ "${typ[$b2]}" = "2" ]
  then
  k=0.8000     # type 2 is conjugated pyridine N atom (C-N-C)  
  else 
  k=1.0000     # type 1 is conjugated carbon atom and type 3 is conjugated pyrrole N atom (C=N-C)
  fi
  n1=`expr \( $b1 - 1 \) \* $atom + $b2`
  n2=`expr \( $b2 - 1 \) \* $atom + $b1`
  matrix[$n1]=`echo "scale=10;$k*$beta"|bc`
  matrix[$n2]=`echo "scale=10;$k*$beta"|bc`
done


for ((i=1;i<=${atom};i++))
do
if [ "${typ[$i]}" != "0" ]
then
  for ((j=1;j<=${atom};j++))
  do
  e=`expr \( $i - 1 \) \* $atom + $j`  
  if [ ! "${matrix[$e]}" ] 
  then
  matrix[$e]=0 
  fi 

  if [ "${typ[$j]}" != "0" ]
  then
  printf "%6.2f" ${matrix[$e]} 
  printf "%16.10f" ${matrix[$e]} >> ./secular_determinant.txt
  fi
  done
printf "\n"
printf "\n" >> ./secular_determinant.txt
fi
done





echo "=============================================================================="
echo "                      CALCULATING ORBITAL ENERGY                              "
echo "=============================================================================="

#Diagonalize secular determinant

echo "import numpy as np                              " >> eigval.py
echo "from numpy import loadtxt                       " >> eigval.py
echo "from numpy import real                          " >> eigval.py
echo "a = loadtxt('secular_determinant.txt')          " >> eigval.py
echo "b = np.linalg.eigvals(a)                        " >> eigval.py
echo "b = real(b)                                     " >> eigval.py
echo "c = np.sort(b)                                  " >> eigval.py
echo "s = ''                                          " >> eigval.py
echo "f = open('eigenvalues.txt','w')                 " >> eigval.py
echo "for line in c:                                  " >> eigval.py
echo "    s = s + str(line) + '\n'                    " >> eigval.py
echo "f.write(s)                                      " >> eigval.py
python ./eigval.py

energy=`awk '{printf "%.10f\n", $1}' eigenvalues.txt`

rm secular_determinant.txt eigenvalues.txt eigval.py

base=`expr ${ele_num} / 2`
resi=`expr ${ele_num} % 2`

if [ "$resi" = "1" ]
then
base=`expr $base + 1`
fi
homo=$base
lumo=`expr $base + 1`

echo "=============================================================================="
printf "%26s%29s\n" "Orbital" "Energy (a.u.)"
for ((a=1;a<=$orb;a++))
do
orb_energy=`echo "$energy"|sed -n "$a"p` 
if [ "$a" = "$homo" ]
then 
orb_type=HOMO
homo_energy=$orb_energy
elif [ "$a" = "$lumo" ]
then 
orb_type=LUMO
lumo_energy=$orb_energy
gap_energy=$(printf "%.10f" `echo "scale=10;${lumo_energy}- ${homo_energy}"|bc`)
degenerated=`echo "$gap_energy<0.000016"|bc`
  if [ "$degenerated" = "1" ]
  then
  lumo=`expr $homo + 1`
  orb_type=HOMO
  fi
else
orb_type=""
fi
printf "%23s%31s%5s\n" "$a" "$orb_energy" "$orb_type"
done

echo "=============================================================================="
echo "=============================================================================="
echo "             HOMO-LUMO ENERGY GAP:  $gap_energy  a.u.                         "
echo "=============================================================================="
echo "=============================================================================="
echo "                       GENERATING SPECTRA DATA                                "
echo "=============================================================================="
rm ./spectra.txt >/dev/null 2>&1
echo "import math                 " >> ./plot.py
echo "homolumogap = $gap_energy   " >> ./plot.py
echo "f = open('spectra.txt','w') " >> ./plot.py
echo "f.write('Wavelength(nm) Intensity\n') " >> ./plot.py
echo "l = ''                      " >> ./plot.py
echo "for x in range(100,1001,2): " >> ./plot.py
echo "    y = math.sqrt(2/math.pi)/70*math.e**(-2*((x-(2.998*10**8)*10**9/((1.60217656*10**(-19))/(6.626*10**(-34)))/(homolumogap))/70)**2) " >> ./plot.py
echo "    l = l + str(x) + '\t '+ str(y) + '\n' "  >> ./plot.py
echo "f.write(l)                  " >> ./plot.py
python ./plot.py
rm ./plot.py
echo "=============================================================================="
echo "                                COMPLETE                                      "
echo "=============================================================================="


