for i in SF/L8Tx*Ty*;do
  Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
  Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
  ./Binder1.x 8 $Tx $Ty $i;
done  > Binder/L8 &

for i in SF/L4Tx*Ty*;do
  Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
  Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
  ./Binder1.x 4 $Tx $Ty $i;
done  > Binder/L4 &

for i in SF/L16Tx*Ty*;do
  Ty=`echo $i | sed 's/Ty/\t/ ' | sed 's/MCS/\t/ ' | cut -f2`;
  Tx=`echo $i | sed 's/Tx/\t/ ' | sed 's/Ty/\t/ ' | cut -f2`;
  ./Binder1.x 16 $Tx $Ty $i;
done > Binder/L16 &
