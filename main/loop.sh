#!/bin/bash

for i in {1..200} 
do

   if [ $i -lt 10 ]; then
   export OUTROOTFILE=outputwithshielding/muons_00$i.root
   elif [ $i -lt 100 ]; then
   export OUTROOTFILE=outputwithshielding/muons_0$i.root
   elif [ $i -lt 1000 ]; then
   export OUTROOTFILE=outputwithshielding/muons_$i.root
   fi

   echo "################################################"
   echo "running simulation $i ..."
   echo "outputwithshielding: $OUTROOTFILE"
   echo "################################################"   
   echo ""

   ./build/b02 runMany.mac

   echo ""
   echo "################################################"
   echo "simulation $i finished"
   echo "################################################"
   echo ""

   
done
