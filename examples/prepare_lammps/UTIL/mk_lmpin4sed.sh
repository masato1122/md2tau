#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

# ------ parameters -----------
#time_npt=2          # ps
#time_pre=2          # ps
#nstructures=260       # 2^8
#nstructures=32770    # 2^15
#nstructures=65540    # 2^16
# -----------------------------
NTHER=10; NDUMP=2

if [ $# -ge 3 ]; then
    time_npt=$1
    time_pre=$2
    nstructures=$3
else
    echo -n " Time for NPT cal [ps]: "
    read time_npt
    echo -n " Time for pre-cal [ps]: "
    read time_pre
    echo -n " # of obtained structures (260, 32770, 65540, etc.): "
    read nstructures
fi

echo ""
echo " # of structures : $nstructures"
echo ""

tmp=`echo "$time_npt / 0.0005" | bc -l`
nptrun=${tmp%.*}

rand=`expr $(( $RANDOM % 100 )) \* $NTHER`
tmp=`echo "$rand + $time_pre / 0.0005" | bc -l`
prerun=${tmp%.*}
if [ $(( $prerun % $NDUMP )) != 0 ]; then
    echo "Error"
    exit
fi

echo ""
echo " NPT : $nptrun steps ( $time_npt ps )"
echo " PRE : $prerun steps ( $time_pre ps + alpha )"
echo ""

for mol in swnt pea; do
OFILE=band_${mol}.in
cat >$OFILE<<EOF
clear

#----------- Parameters --------------------------
variable tstep     equal 0.0005
variable tdamp     equal 0.05
variable TEMP      equal 300.00

EOF
if [ $nptrun != "0" ]; then
cat >>$OFILE<<EOF
variable NPTRUN    equal ${nptrun}
variable NPTTHER   equal \${NPTRUN}/${NTHER}
variable NPTDUMP   equal \${NPTRUN}/${NDUMP}

EOF
fi
cat >>$OFILE<<EOF
variable PRERUN    equal ${prerun}
variable PRETHER   equal \${PRERUN}/${NTHER}
variable PREDUMP   equal \${PRERUN}/${NDUMP}

variable DUMPSTEP  equal 16
variable NRUN      equal ${nstructures}*\${DUMPSTEP}
variable THERSTEP  equal \${NRUN}/10

variable kb        equal 1.3806488e-23
variable ee        equal 1.60217657e-19

# ---------- Parameters for LAMMPS ---------------
units         metal
atom_style    atomic
dimension     3
boundary      p p p
atom_modify   map array
read_data     data.lammps
timestep      \${tstep}

EOF

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ---- for SWNT or peapod
if [ $mol == "swnt" ]; then
cat >>$OFILE<<EOF
# ---------- Define Interatomic Potential --------
pair_style tersoff
pair_coeff * * opt.tersoff C
EOF
else
cat >>$OFILE<<EOF
# ---------- Define Interatomic Potential --------
variable LJcutoff equal 12.00
variable epsilon  equal 0.00239
variable sigma    equal 5.0

pair_style hybrid lj/cut \${LJcutoff} tersoff tersoff tersoff
pair_coeff * * tersoff 1 opt.tersoff C NULL NULL
pair_coeff * * tersoff 2 opt.tersoff NULL C NULL
pair_coeff * * tersoff 3 opt.tersoff NULL NULL C
pair_coeff 1 2 lj/cut \${epsilon} \${sigma}
pair_coeff 1 3 lj/cut \${epsilon} \${sigma}
pair_coeff 2 3 lj/cut \${epsilon} \${sigma}
EOF
fi
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cat >>$OFILE<<EOF

neighbor      0.3 bin
neigh_modify  delay 0

#--------- Start the calculation -------------

reset_timestep 0

EOF

# >>>>>
# --------- if nptrun==0, minimize the structure
if [ $nptrun == "0" ]; then
cat >>$OFILE<<EOF
#--- 1. minimization
fix         1 all box/relax z 10.0
dump        1 all custom 100000 mini.dump id type x y z
dump_modify 1 sort id
dump_modify 1 format float %.5f
thermo      1000
minimize    1e-20 0.00001 50000 50000
unfix       1
undump      1

EOF

# --------- if nptrun==0, perform NPT 
else
cat >>$OFILE<<EOF
#--- 1. NPT
velocity all create 300.0 ${RANDOM} rot yes dist gaussian
fix      1 all npt temp 300 300 \${tdamp} z 0.0 0.0 0.1
thermo   \${NPTTHER}

dump        id1 all custom \${NPTDUMP} npt.dump id type x y z
dump_modify id1 sort id
dump_modify id1 format float %.5f

run      \${NPTRUN}
unfix    1
undump   id1

EOF
fi
# <<<<<

cat >>$OFILE<<EOF
#--- 2. pre-MD
reset_timestep 0

velocity  all create \${TEMP} ${RANDOM} rot yes dist gaussian
fix       2 all nvt temp 300. 300. \${tdamp}

thermo    \${PRETHER}
dump        2 all custom \${PREDUMP} pre.dump id type x y z vx vy vz
dump_modify 2 sort id
dump_modify 2 format float %.5f
run         \${PRERUN}
unfix       2
undump      2

#--- 3. MD for SED
reset_timestep 0

fix         3 all nvt temp 300. 300. \${tdamp}
thermo      \${THERSTEP}
dump        3 all custom \${DUMPSTEP} md4band.dump id type x y z vx vy vz
dump_modify 3 sort id
dump_modify 3 format float %.5e
run         \${NRUN}
unfix       3
undump      3

print "All done!"
EOF
done

