#!/bin/bash

cd $1
mkdir layer_comp
for i in layer0/*; do
  picname=${i#layer0/}
  echo 'processing ' $picname
  convert -append -gravity center -background LightGray -size x40 layer0/$picname label:'1025.6' layer1/$picname label:'1026.85' layer2/$picname label:'1027.4' layer3/$picname label:'1027.7' layer4/$picname label:'1028.1' layer_comp/$picname

done

cd ..

