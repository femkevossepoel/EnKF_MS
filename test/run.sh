export OMP_NUM_THREADS=14


for i in infile.*
do
   echo $i
   sep=$(echo ${i:(-3)})
   exp=$(echo "${i/infile./}")
   expinfile=$(sed -n '2p' $i | sed -e 's/  */ /g' | cut -f2 -d' ')
   [[ $exp != $expinfile ]] && sed -i "s/$expinfile/$exp/g" $i

   if [ -f $i ]
   then
      cp $i infile.in
      multiscale   | tee  ${exp}.log
      RESULT=$?
      if [ $RESULT == 0 ]; then
         pushd $exp
         gnuplot ../cpdf2.gnu
         gnuplot ../rms.gnu
         [[ -f costf.dat ]] &&  gnuplot ../costf.gnu
         popd
      else
        echo failed 2
      fi
   fi
done
