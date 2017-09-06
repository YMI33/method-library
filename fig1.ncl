;************************************************
; Basic User input
;************************************************

diri = "/data1/JRA55/monthly/pressure/"
time_s=277 ;1981/1 for JRA55
N=30 ; total count of year
mon=8
plot_format="X11"
plot_name="fig1"
plot_color="BlueDarkRed18"          ;"MPL_YlOrBr"       ;"BlueWhiteOrangeRed"

;************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./weighted_composite_analysis.ncl"
  

begin
;-------------------------read data-----------------------------
;index
fileindex=addfile("./augQ1.nc","r")
Q1=doubletofloat(fileindex->Q1)
Q18=Q1(0:29)

;field
fili01="uwnd.jra55.mon.mean.nc"
fili02="vwnd.jra55.mon.mean.nc"
fili03="hgt.jra55.mon.mean.nc"
fin01=addfile(diri+fili01,"r")
fin02=addfile(diri+fili02,"r")
fin03=addfile(diri+fili03,"r")

var01=fin01->uwnd(time_s-1:time_s+N*12-2,{850},{-20:40},{90:180})
var02=fin02->vwnd(time_s-1:time_s+N*12-2,{850},{-20:40},{90:180})
var03=fin03->hgt(time_s-1:time_s+N*12-2,{500},{-20:40},{90:180})

;time & geographic information
lat03=var03&lat
lon03=var03&lon
dim_var=dimsizes(var03)
var03:=reshape(var03,(/N,12,dim_var(1),dim_var(2)/))
var03:=var03(:,mon-1,:,:)
var02:=reshape(var02,(/N,12,dim_var(1),dim_var(2)/))
var02:=var02(:,mon-1,:,:)
var01:=reshape(var01,(/N,12,dim_var(1),dim_var(2)/))
var01:=var01(:,mon-1,:,:)

var03!0="year"
var03!1="lat"
var03!2="lon"
var03&lon=lon03
var03&lat=lat03
copy_VarCoords(var03,var01)
copy_VarCoords(var03,var02)

;--------------------------------WCA--------------------------------
var_wca=var03(lat|:,lon|:,year|:)
output=weighted_composite_analysis(Q18,var_wca({-20:40},{90:180},:),200,0.95)
wca_hgt=output[0]
mc_test_hgt=output[1]

var_wca=var01(lat|:,lon|:,year|:)
output=weighted_composite_analysis(Q18,var_wca({-20:40},{90:180},:),200,0.95)
wca_uwnd=output[0]
mc_test_uwnd=output[1]

var_wca=var02(lat|:,lon|:,year|:)
output=weighted_composite_analysis(Q18,var_wca({-20:40},{90:180},:),200,0.95)
wca_vwnd=output[0]
mc_test_vwnd=output[1]

;----------------------------plot-------------------------------------
wca_test_f_hgt=where(mc_test_hgt,wca_hgt,wca_hgt@_FillValue)
copy_VarMeta(wca_hgt,wca_test_f_hgt)
wca_test_f_uwnd=where(mc_test_uwnd.or.mc_test_vwnd,wca_uwnd,wca_uwnd@_FillValue)
copy_VarMeta(wca_uwnd,wca_test_f_uwnd)
wca_test_f_vwnd=where(mc_test_uwnd.or.mc_test_vwnd,wca_vwnd,wca_vwnd@_FillValue)
copy_VarMeta(wca_vwnd,wca_test_f_vwnd)

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
res@cnMaxlevelValF=20
res@cnMinLevelValF=-20

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

;resWind
resWind=True
resWind@gsnFrame=False
resWind@gsnDraw=False
resWind@gsnSpreadColors      = True

resWind@lbLabelBarOn=False
resWind@vcRefMagnitudeF=1
resWind@vcRefLengthF=0.01
resWind@vcMinDistanceF               =0.03    ; density of vector 
;resWind@vcGlyphStyle="CurlyVector"
;resWind@vcLevelSelectionMode= "ExplicitLevels"
;resWind@vcLevels= (/90/)
resWind@vcLevelColors= (/ "grey","black"/)
;resWind@vcMinDistanceF=0.015
opt = True
opt@gsnShadeFillType = "pattern"      ; pattern fill
opt@gsnShadeHigh = 17                   ; use pattern #2 

indicate=2

plot=gsn_csm_contour_map(wks,wca_test_f_hgt(:,:,indicate),res)
plot_shade=gsn_csm_contour(wks,wca_test_f_hgt(:,:,indicate),resShade)
plot_test=gsn_contour_shade(plot_shade,-5,5,opt)
overlay(plot,plot_test)
plot_vec=gsn_csm_vector(wks,wca_test_f_uwnd(:,:,indicate),wca_test_f_vwnd(:,:,indicate),resWind)
overlay(plot,plot_vec)

draw(plot)
frame(wks)


end
