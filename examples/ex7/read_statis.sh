a=0
while read line
do
        a=$((a+1))
        if [ $a -eq "9" ];then
                a=0
                echo $step $vdw
        fi
        if [ $a -eq "1" ];then
                step=`echo $line | awk '{print $1}'`
        fi
        if [ $a -eq "2" ];then
                vdw=`echo $line | awk '{print $4}'`
        fi
#       echo $a $line   
done < STATIS

