for v in $(seq 0.57 0.01 0.63); do
  for i in SF/L*Tx10000Ty*;do
    Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
    Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
    L=`echo $i | sed 's/L/\t/ ' | sed 's/Tx/\t/ ' | cut -f2`;
    ./Fine1.x $L $Tx $Ty 3.2 $v $i;
  done  > Fine1/1xx_2xy_3s10_4s01_5y_6L_710000_8Ty_9Tc_10v
done &

for v in $(seq 0.57 0.01 0.63); do
  for i in SF/L*Tx5Ty*;do
    Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
    Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
    L=`echo $i | sed 's/L/\t/ ' | sed 's/Tx/\t/ ' | cut -f2`;
    ./Fine1.x $L $Tx $Ty 3.2 $v $i;
  done  > Fine1/1xx_2xy_3s10_4s01_5y_6L_75_8Ty_9Tc_10v
done &

for v in $(seq 0.57 0.01 0.63); do
  for i in SF/L*Tx3.5Ty*;do
    Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
    Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
    L=`echo $i | sed 's/L/\t/ ' | sed 's/Tx/\t/ ' | cut -f2`;
    ./Fine1.x $L $Tx $Ty 3.2 $v $i;
  done  > Fine1/1xx_2xy_3s10_4s01_5y_6L_73.5_8Ty_9Tc_10v
done &