#PBS -S /bin/bash
#PBS -q batch
#PBS -N prepare_DMR
#PBS -l nodes=1:ppn=2:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=30gb
#PBS -j oe

cd $PBS_O_WORKDIR

#Setup DMR files
mkdir DMR
mkdir DMR/allc DMR/CG DMR/CH DMR/CN
for i in FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3
do
  echo "Copying" "$i"
  cp "$i"/methylCseq/allc/"$i"_allc.tar.bz2 DMR/allc/
done

cd DMR/allc
for i in *tar.bz2
do
  echo "Unpacking" "$i"
  tar -xjvf "$i"
done
cd ../../
