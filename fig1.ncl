;************************************************
; Basic User input
;************************************************

diri = "/data1/JRA55/monthly/pressure/"
time_s=277 ;1981/1 for JRA55
N=33 ; total count of year
mon=8
plot_format="X11"
plot_name="fig1"
plot_color="BlueDarkRed18"          ;"MPL_YlOrBr"       ;"BlueWhiteOrangeRed"

;************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin
;-------------------------read data-----------------------------
;fileindex=addfile("./Q1_eof_aug_1981-2013.nc","r")
;Q18=fileindex->eof_ts(0,:)
fileindex=addfile("./augQ1.nc","r")
Q18=doubletofloat(fileindex->Q1)

fili03="hgt.jra55.mon.mean.nc"


fin03=addfile(diri+fili03,"r")
;printVarSummary(fin03->hgt)
var03=fin03->hgt(time_s-1:time_s+N*12-2,{500},{-20:40},{90:180})
lat03=var03&lat
lon03=var03&lon
dim_var=dimsizes(var03)
var03:=reshape(var03,(/N,12,dim_var(1),dim_var(2)/))
var03:=var03(:,mon-1,:,:)

var03!0="year"
var03!1="lat"
var03!2="lon"
var03&lon=lon03
var03&lat=lat03

var_wca=var03(lat|:,lon|:,year|:)

;===================================================================
;----------------------------WCA-------------------------------------
;===================================================================
samplesize=200
alpha=0.95
windex=Q18-avg(Q18)

;dimension information
dimsiz=dimsizes(var_wca)
siznum=dimsizes(dimsiz)
wca=new(array_append_record(dimsiz(0:1),(/3/),0),typeof(var_wca),var_wca@_FillValue)

;calculate the weight
ind_p=ind(windex.gt.0)
ind_n=ind(windex.lt.0)
N=dimsiz(siznum-1)
totalweight=sum(windex(ind_p))
w_multidim=conform(var_wca,windex,siznum-1)

;weighted composite analysis
wca(:,:,0)=dim_sum_Wrap(var_wca(:,:,ind_p)*w_multidim(:,:,ind_p))/totalweight
wca(:,:,1)=-dim_sum_Wrap(var_wca(:,:,ind_n)*w_multidim(:,:,ind_n))/totalweight
wca(:,:,2)=wca(:,:,0)-wca(:,:,1)


copy_VarMeta(var_wca,wca)
wca!2="mode"
wca&mode=(/"positive","negative","difference"/)

print("start test progress")

;preset for MC test
wca_set=new(array_append_record(dimsiz(0:1),(/3,samplesize/),0),typeof(var_wca),var_wca@_FillValue)
rseed1 = toint(systemfunc(" date +%s"))%100000
rseed2 = toint(103*rseed1)

do i=0,samplesize-1

; generate new order
rseed1=rseed1+i
rseed2=rseed2+i
random_setallseed(rseed1, rseed2)
neworder=generate_unique_indices(N)
ind_new=windex(neworder)

;regenerate procedure
ind_p=ind(ind_new.gt.0)
ind_n=ind(ind_new.lt.0)
w_multidim=conform(var_wca,ind_new,siznum-1)

wca_set(:,:,0,i)=dim_sum_Wrap(var_wca(:,:,ind_p)*w_multidim(:,:,ind_p))/totalweight
wca_set(:,:,1,i)=-dim_sum_Wrap(var_wca(:,:,ind_n)*w_multidim(:,:,ind_n))/totalweight
wca_set(:,:,2,i)=wca_set(:,:,0,i)-wca_set(:,:,1,i)

end do

;-------------------obtain confidence interval-------------------------

boun_left=floattoint((1-alpha)*samplesize)-1
boun_right=floattoint(alpha*samplesize)-1
mc_test_array=wca_set
mc_test_array_ind=dim_pqsort(mc_test_array,2)

;-----------------------test WCA result---------------------------------
mc_test=(wca.le.mc_test_array(:,:,:,boun_left)).or.(wca.ge.mc_test_array(:,:,:,boun_right))
;========================== WCA END ===================================
;======================================================================

;----------------------------plot-------------------------------------
wca_test_f=where(mc_test,wca,wca@_FillValue)
copy_VarMeta(wca,wca_test_f)

wks=gsn_open_wks(plot_format,plot_name)     
gsn_define_colormap(wks,plot_color)

res=True
res@gsnFrame=False
res@gsnDraw=False
res@gsnSpreadColors      = True
res@cnFillOn=True
res@cnLinesOn=False
res@cnLineLabelsOn=False
resShade=res
res@gsnAddCyclic=False
res@mpMinLatF            = -20     ;24   ; zoom in on map
res@mpMaxLatF            = 40
res@mpMinLonF            = 90
res@mpMaxLonF            = 180

res@vpWidthF=0.8
res@vpHeightF=0.5
res@gsnMaximize=True
res@tmXBMode="Explicit"
;res@tmXBValues=ispan(110,130,1)
res@lbLabelBarOn=True

;resShade
resShade@cnFillOn=False
resShade@cnLinesOn=False
resShade@cnLineLabelsOn=False
resShade@cnInfoLabelOn=False

opt = True
opt@gsnShadeFillType = "pattern"      ; pattern fill
opt@gsnShadeHigh = 17                   ; use pattern #2 


plot=gsn_csm_contour_map(wks,wca(:,:,2),res)
plot_shade=gsn_csm_contour(wks,wca_test_f(:,:,2),resShade)
plot_test=gsn_contour_shade(plot_shade,-5,5,opt)
overlay(plot,plot_test)
draw(plot)
frame(wks)
;plot=gsn_csm_contour_map(wks,wca_test(:,:,0),res)




end

