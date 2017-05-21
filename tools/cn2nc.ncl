
; Description: convert rgm output txt data to nc readable by grads
;
;      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
;     Created: 2017-05-17 09:18:04 BJT
; Last Change: 2017-05-17 10:57:46 BJT

;in h15.txt:
;txt, nrow = 140 (lat), ncol = 140 (lon),

; data info <<<1

data_dir = "/home/ouniansen/data/rgm"
inputfile = data_dir + "/h15.txt"
fileOut = data_dir + "/h15.nc"
nlat = 140
nlon = 140
nrec = 1

; retrieve data from .dat <<<1

h = new((/nrec,nlat,nlon/), float)

do rec = 0, nrec - 1
    h(rec,:,:) = asciiread(inputfile,(/nlat,nlon/),"float")
end do

; assign atts <<<1

h = h*0.01
h@units = "m"
h@long_name    = "layer thickness"

h!0 = "time"
h!1 = "lat"
h!2 = "lon"

lat = fspan(-69,70,nlat) ; NOT ispan
lat!0 = "lat"
lat@units = "degrees_north"

lon = fspan(0,139,nlon)
lon!0 = "lon"
lon@units = "degrees_east"

time = ispan(15,15,1)
time@long_name = "time"
time@units = "years since 0000-01-01"

h&time = time
h&lat = lat
h&lon = lon

system("rm -f " + fileOut)
f = addfile(fileOut, "c")
f->h = h
f@source = "create by " + systemfunc("pwd") + "/" + getenv("NCL_ARG_0") + \
", " + systemfunc("date")

delete(f)
