#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
if [ $# -ge 1 ]; then
    ical=$1
else
    echo -n "Input calculation ID: "
    read ical
fi

WORKSEKI=/work/k0437/k043700/Peapod

OFILE1=1_fermi.sh
OFILE2=2_seki.sh
OFILE3=3_fermi.sh

# ---- 1. for fermi
nprocs=20
cat >$OFILE1<<EOF
#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=${nprocs}
#PBS -j oe
#PBS -N CAL${ical}_md

export LANG=C
export OMP_NUM_THREADS=1
cd \$PBS_O_WORKDIR
nprocs=${nprocs}
MPIRUN=mpirun

cd \$PBS_O_WORKDIR
for mol in swnt pea; do
    DIR0=\$mol
    DIR1=\${mol}.tar.gz
    if [ ! -e \$DIR0 ]; then
        echo "Error \$DIR0 does not exist."
        exit
    fi
    
    cd ./\$DIR0
    \$MPIRUN -np \${nprocs} lmp_intel < band_\${mol}.in > log_md.txt
    
    cd ../
    tar -zcvf \${DIR1} \${DIR0}
done

# ---- get name of directory
line=\`pwd\`
data=(\`echo \$line | tr -s '/' ' '\`)
n=\${#data[*]}; i=\`expr \$n - 1\`
dname=\${data[\$i]}

# --- send directory to sekirei
cd \$PBS_O_WORKDIR
sekirei="sekirei.issp.u-tokyo.ac.jp"
DGO=${WORKSEKI}/\${dname}
ssh -l k043700 \${sekirei} mkdir -p \$DGO
for mol in swnt pea; do
    scpsekirei_send2 \${mol}.tar.gz \$DGO
done
scpsekirei_send2 2_seki.sh \$DGO

# --- run
ssh -l k043700 \${sekirei} qsub \$DGO/${OFILE2}

# --- remove original directories
cd \$PBS_O_WORKDIR
for mol in swnt pea; do
    rm -r \${mol}
done
EOF

# --- for sekirei
WORK_seki=/work/k0437/k043700
cat >$OFILE2<<EOF
#!/bin/sh
#QSUB -queue F2fat
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=20:00:00
#PBS -N CAL${ical}_sed

cd \$PBS_O_WORKDIR
nprocs=24

# ---- may need to be modified
DSEK=${WORKSEKI}/CAL${ical}

if [ ! -e \$DSEK ]; then
    echo "Error: cannot find \$DSEK"
    exit
fi

LOG=log_sed.txt
cd \$DSEK
for mol in swnt pea; do
    DIR0=\${mol}
    DIR1=\${mol}.tar.gz
    DIR2=\${mol}_done.tar.gz
    
    cd \$DSEK

    # --- expand & remove
    tar -zxvf \$DIR1
    
    # --- calculate SED
    cd ./\$DIR0
    date > \$LOG
    mpirun -np 5 dump2sed sed_all.in >> \$LOG
    date >> \$LOG
    rm md4band.dump
    
    # --- compression
    cd ../
    tar -zcvf \${DIR2} \${DIR0}
    rm \${DIR1}
done

# ---- get name of directory and send directory
#for mol in swnt pea; do
#    scpfermi_send2 ./\${mol}_done.tar.gz \$DFER
#done
# --- run the job in Fermi
#ssh -p 8656 -l ohnishi fermi.t.u-tokyo.ac.jp qsub \$DFER/3_fermi.sh
EOF

# --- 3. get result from sekirei and expand
cat >$OFILE3<<EOF
#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=2
#PBS -j oe
#PBS -N CAL${ical}_expand

export LANG=C
export OMP_NUM_THREADS=1
cd \$PBS_O_WORKDIR
nprocs=2

cd \$PBS_O_WORKDIR

# ---- get name of directory
line=\`pwd\`
data=(\`echo \$line | tr -s '/' ' '\`)
n=\${#data[*]}; i=\`expr \$n - 1\`
dname=\${data[\$i]}

for mol in swnt pea; do
    DIR0=\${mol}
    DIR1=\${mol}.tar.gz
    DIR2=\${mol}_done.tar.gz
    for d in \$DIR0 \$DIR1; do
        if [ -e \$d ]; then
            mv \$d X_\${d}
        fi
    done
    scpsekirei_get ${WORKSEKI}/\${dname}/\${DIR2}
    tar -zxvf \${DIR2}
done

# --- change the directory name
cd ../
mv ./\${dname} ./\${dname}_O
EOF

