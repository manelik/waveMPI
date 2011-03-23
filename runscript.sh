

#Usage
#
# runscript.sh numproc binary input

echo 'Script '$0' running'
mpirun -n $1 $2 < $3