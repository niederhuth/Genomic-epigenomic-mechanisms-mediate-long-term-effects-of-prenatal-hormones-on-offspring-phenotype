#PBS -S /bin/bash
#PBS -q batch
#PBS -N browser_files
#PBS -l nodes=1:ppn=10:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=20gb
#PBS -j oe

echo "Starting"
cd $PBS_O_WORKDIR
module load python/3.5.1
mkdir browser_files

for i in FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 \
MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3
do
  echo "Making browser files for" "$i"
  cd "$i"/methylCseq/allc
  tar -xjvf "$i"_allc_total.tar.bz2
  python3.5 /home/chadn/bin/allc2bigwig/allc_to_bigwig.py -all -sort -p=10 \
  ../../../ref/bowtie2/Tguttata_v3.2.4.fa.fai "$i"_allc_total.tsv
  cd ../../../
done
cp */methylCseq/allc/*.bw.c* browser_files/

echo "Done"
