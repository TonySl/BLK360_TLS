# Shengli Tao
# Toulouse, France
# sltao1990@gmail.com

import os,glob

indir=r'E:\TLidar\2019\StripALL_laz_optimize_thinCC1cm_clipboundary'

os.chdir(indir)
filenames = glob.glob("*.laz")

f = open(r"E:\TLidar\2019\CC_subsample_1cm_singlescan_addsourceID.bat","w")


for i in range(0,len(filenames)):
        print filenames[i]
        tif1= filenames[i]
##      CloudCompare -O E:\TLidar\2019\StripALL_laz_optimize\StripALL_0_1.laz -C_EXPORT_FMT LAS -SS SPATIAL 0.01  -CLEAR
        line2="las2las -i E:\\TLidar\\2019\\StripALL_laz_optimize_thinCC1cm_clipboundary\\" + tif1 + " " + "-o"+ " " + "E:\\TLidar\\2019\\StripALL_laz_optimize_thinCC1cm_clipboundary_withscanID\\"+tif1 +" "+"-change_point_source_from_to 1"+" "+ str(i+1)
        
        f.write(line2+"\n")

       
f.close()

##os.system("CC_subsample_1cm_singlescan_addsourceID.bat")
