

#Usage
#
# runscript.sh numproc binary input

if [ ${#@} -lt 3 ]
    then
    echo 'Bad usage'
    echo $0' numproc binary input'
    exit 0
fi

echo 'Script '$0' running'
mpirun -n $1 $2 < $3