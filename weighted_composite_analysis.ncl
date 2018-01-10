;***********************************************************************
;NCL script for weighted composite analysis (WCA), written by Xie Zhiang, 2017/9/6
;A NCL function for WCA with Mote Carlo two-tail test
;You can use the function by output=weighted_composite_analysis(index,field,N,alpha)
:where index is an index used as standard to tell positive and negative event, 
;field is a field ready for composition, which must be 3D variable and the rightmost dimension must be same as length of index
;N is the ensemble size of Mote Carlo test, 200 and above is recommended
;alpha is the significant level for test
;output contains two variable:
;output[0] is the result of WCA, a 3D variable and you can find positive (0), negative (1), and difference (2) in rightmost demension, 
;output[1] is the result of Monte Carlo two-tail test, +1 & -1 represent passing the significant test, 0 represents not passing the test. 

;***********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("weighted_composite_analysis")
function weighted_composite_analysis(index_input,var_wca,samplesize,alpha0)
local wca,mc_test


begin
;===================================================================
;----------------------------WCA-------------------------------------
;===================================================================


windex=index_input-avg(index_input)

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

print("test progress start")

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
alpha=(1+alpha0)/2
boun_left=floattoint((1-alpha)*samplesize)-1
boun_right=floattoint(alpha*samplesize)-1
mc_test_array=wca_set
mc_test_array_ind=dim_pqsort(mc_test_array,2)

;-----------------------test WCA result---------------------------------
mc_test=(wca.le.mc_test_array(:,:,:,boun_left)).or.(wca.ge.mc_test_array(:,:,:,boun_right))

print("test progress end")

return([/wca,mc_test/])

end
