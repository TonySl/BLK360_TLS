# Shengli Tao
# Toulouse, France
# sltao1990@gmail.com

import os,glob

indir=r'E:\TLidar\2019\StripALL_laz_optimize'

os.chdir(indir)
filenames = glob.glob("*.laz")

f = open(r"E:\TLidar\2019\CC_subsample_1cm_singlescan.bat","w")


for i in range(0,len(filenames)):
        print filenames[i]
        tif1= filenames[i]
##      CloudCompare -O E:\TLidar\2019\StripALL_laz_optimize\StripALL_0_1.laz -C_EXPORT_FMT LAS -SS SPATIAL 0.02  -CLEAR
        line2="CloudCompare -SILENT -O E:\\TLidar\\2019\\StripALL_laz_optimize\\" + tif1 + " " + "-C_EXPORT_FMT" + " " + "LAS" + " " + "-SS" + " " + "SPATIAL" + " " + "0.02" + " " + "-CLEAR"
        
        f.write(line2+"\n")

       
f.close()

##os.system("CC_subsample_1cm_singlescan.bat")
