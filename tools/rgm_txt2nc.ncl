
; Description: convert rgm output txt data to nc readable by grads
;
;      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
;     Created: 2017-05-17 09:18:04 BJT
; Last Change: 2017-05-17 22:42:29 BJT

;txt, nrow = 140 (lat), ncol = 140 (lon),
; 15 + 1 times, one time per year.

; data info <<<1

data_dir = "/home/ouniansen/rgm/Pacific"
inputfile = data_dir + "/fort.203"
fileOut = data_dir + "/post/v.nc"
nlat = 140
nlon = 140
nrec = 16

; retrieve data <<<1

v = new((/nrec,nlat,nlon/), float)
data = new((/nlat*nrec,nlon/), float)

data = asciiread(inputfile,(/nlat*nrec,nlon/),"float")
do rec = 1, nrec
    row1 = (rec-1)*nlat
    row2 = row1 + nlat - 1
    v(nrec-1,:,:) = data(row1:row2,:)
end do

; assign atts <<<1

v = v*0.01
v@units = "m/s"
v@long_name = "zonal velocity"

v!0 = "time"
v!1 = "lat"
v!2 = "lon"

lat = fspan(-69,70,nlat) ; NOT ispan
lat!0 = "lat"
lat@units = "degrees_north"

lon = fspan(0,139,nlon)
lon!0 = "lon"
lon@units = "degrees_east"

time = ispan(0,15,1)
time@long_name = "time"
time@units = "years since 0000-01-01"

v&time = time
v&lat = lat
v&lon = lon

system("rm -f " + fileOut)
f = addfile(fileOut, "c")
f->v = v
f@source = "create by " + systemfunc("pwd") + "/" + getenv("NCL_ARG_0") + \
", " + systemfunc("date")

delete(f)
