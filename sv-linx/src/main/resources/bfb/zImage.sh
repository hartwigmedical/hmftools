ImageMagick

Trim
convert -density 300 cell.pdf[1]  -flatten -fuzz 10% -trim +repage allphase.png
convert -density 300 cell.pdf[2]  -flatten -fuzz 10% -trim +repage m1.png
convert -density 300 cell.pdf[3]  -flatten -fuzz 10% -trim +repage m2.png
convert -density 300 cell.pdf[4]  -flatten -fuzz 10% -trim +repage g2.png
convert -density 300 cell.pdf[5]  -flatten -fuzz 10% -trim +repage g1.png
convert -density 300 cell.pdf[6]  -flatten -fuzz 10% -trim +repage s.png

Overlap Image
convert COLO829T.brb.step00.png ./cell/allphase.png  -geometry +100+100 -composite -density 300  out00.png
convert COLO829T.brb.step01.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out01.png
convert COLO829T.brb.step02.png ./cell/m2.png  -geometry +100+100 -composite -density 300  out02.png
convert COLO829T.brb.step03.png ./cell/s.png  -geometry +100+100 -composite -density 300  out03.png
convert COLO829T.brb.step04.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out04.png

convert COLO829T.brb.step05.png ./cell/m1.png  -geometry +100+100 -composite -density 300  out05.png
convert COLO829T.brb.step06.png ./cell/m2.png  -geometry +100+100 -composite -density 300  out06.png
convert COLO829T.brb.step07.png ./cell/s.png  -geometry +100+100 -composite -density 300  out07.png
convert COLO829T.brb.step08.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out08.png

convert COLO829T.brb.step09.png ./cell/m1.png  -geometry +100+100 -composite -density 300  out09.png
convert COLO829T.brb.step10.png ./cell/m2.png  -geometry +100+100 -composite -density 300  out10.png
convert COLO829T.brb.step11.png ./cell/g1.png  -geometry +100+100 -composite -density 300  out11.png
convert COLO829T.brb.step12.png ./cell/s.png  -geometry +100+100 -composite -density 300  out12.png
convert COLO829T.brb.step13.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out13.png

convert COLO829T.brb.step14.png ./cell/m1.png  -geometry +100+100 -composite -density 300  out14.png
convert COLO829T.brb.step15.png ./cell/m2.png  -geometry +100+100 -composite -density 300  out15.png
convert COLO829T.brb.step16.png ./cell/s.png  -geometry +100+100 -composite -density 300  out16.png
convert COLO829T.brb.step17.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out17.png
convert COLO829T.brb.step18.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out18.png
convert COLO829T.brb.step19.png ./cell/g2.png  -geometry +100+100 -composite -density 300  out19.png

convert COLO829T.brb.step20.png ./cell/m1.png  -geometry +100+100 -composite -density 300  out20.png
convert COLO829T.brb.step21.png ./cell/m2.png  -geometry +100+100 -composite -density 300  out21.png
convert COLO829T.brb.step22.png ./cell/g1.png  -geometry +100+100 -composite -density 300  out22.png
convert COLO829T.brb.step23.png ./cell/allphase.png  -geometry +100+100 -composite -density 300  out23.png



Make GIF
convert -background white -alpha remove -layers OptimizePlus -delay 3x1 out*.png -loop 0  bfb.gif