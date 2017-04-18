#PBS -S /bin/bash
#PBS -q batch
#PBS -N CG_dmr
#PBS -l nodes=1:ppn=22:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=100gb
#PBS -j oe

cd $PBS_O_WORKDIR
echo "Starting"
module load python/2.7.8
module load R/3.1.2
module load gsl/1.16/gcc/4.4.7

#Bdistachyon
echo "Starting B. distachyon DMR analysis"

#CG
cd Bdistachyon/DMR/CG
echo "Calling CG DMRs"
python ../../../../scripts/DMR_scripts/Bdistachyon_CG.py > CG_dmrfind_info.txt

#CHG
cd ../CHG
echo "Calling CHG DMRs"
python ../../../../scripts/DMR_scripts/Bdistachyon_CHG.py > CHG_dmrfind_info.txt

#CHH
cd ../CHH
echo "Calling CHH DMRs"
python ../../../../scripts/DMR_scripts/Bdistachyon_CHH.py > CHH_dmrfind_info.txt

#CNN
cd ../CNN
echo "Calling CNN DMRs"
python ../../../../scripts/DMR_scripts/Bdistachyon_CNN.py > CNN_dmrfind_info.txt
echo "B. distachyon DMR analysis complete"
cd ../../../

#Osativa
echo "Starting O. sativa DMR analysis"

#CG
cd Osativa/DMR/CG
echo "Calling CG DMRs"
python ../../../../scripts/DMR_scripts/Osativa_CG.py > CG_dmrfind_info.txt

#CHG
cd ../CHG
echo "Calling CHG DMRs"
python ../../../../scripts/DMR_scripts/Osativa_CHG.py > CHG_dmrfind_info.txt

#CHH
cd ../CHH
echo "Calling CHH DMRs"
python ../../../../scripts/DMR_scripts/Osativa_CHH.py > CHH_dmrfind_info.txt

#CNN
cd ../CNN
echo "Calling CNN DMRs"
python ../../../../scripts/DMR_scripts/Osativa_CNN.py > CNN_dmrfind_info.txt
echo "O. sativa DMR analysis complete"
cd ../../../
echo "Done"
