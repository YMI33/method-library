function [difference sigcrit]=wcomposite_para(x,y,alpha0)
%计算加权合成分析的函数，采用双侧蒙特卡洛检验,并行版本
%函数形式 [WCA significance]=wcomposite(index,field,alpha0)
%其中index为一维指数，field为加权合成分析的物理量场，必须为三维变量，第三维为时间维，alpha0是统计检验标准，如0.95表示95%统计检验
%输出量WCA为三维变量，前两维与输入的field相同，第三维为不同异常结果，正异常结果为WCA(:,:,1),负异常为WCA(:,:,2)，合成差为WCA(:,:,3)
%significance为显著检验结果，significance(:,:,1)为置信上界，significance(:,:,2)为置信下界,WCA值高于上界或低于下界的即为通过检验
sizlen=size(x);
if sizlen(1)==1
    x=x';
    sizlen=sizlen(2);
else
    sizlen=sizlen(1);
end

%双侧检验
alpha0=(1+alpha0)/2.;

sizy=size(y);
siznum=3;
%计算合成结果
x_anomal=x-mean(x);
y_mean=mean(y,siznum);
y_anomal=y-repmat(y_mean,[ones(1,siznum-1) sizy(siznum)]);
p=find(x_anomal>0);
n=find(x_anomal<0);

%difference第三维分别为：1.positive part; 2.negative part; 3. differnece
difference=nan(sizy(1),sizy(2),3);
for k1=1:sizy(1)
    for k2=1:sizy(2)
      difference(k1,k2,1)=sum(x_anomal(p).*squeeze(y_anomal(k1,k2,p)))/sum(abs(x_anomal(p)))+y_mean(k1,k2);
      difference(k1,k2,2)=-sum(x_anomal(n).*squeeze(y_anomal(k1,k2,n)))/sum(abs(x_anomal(n)))+y_mean(k1,k2);
    end
end
difference(:,:,3)=difference(:,:,1)-difference(:,:,2);

%计算统计检验
N=200;%蒙特卡洛随机试验次数
diff_test=nan(sizy(1),sizy(2),N,3);
%计算随机序列的加权合成分析结果
parfor k=1:N
%   k
  siz_a=size(y_anomal);
  order=randperm(sizlen);
  x_test=x(order);
  x_test=x_test-mean(x_test);
  p=find(x_test>0);
  n=find(x_test<0);
  pos_part=nan(siz_a(1),siz_a(2));
  neg_part=pos_part;
    for k1=1:siz_a(1)
        for k2=1:siz_a(2)
         pos_part(k1,k2)=sum(x_test(p).*squeeze(y_anomal(k1,k2,p)))/sum(abs(x_test(p)));
         neg_part(k1,k2)=-sum(x_test(n).*squeeze(y_anomal(k1,k2,n)))/sum(abs(x_test(n)));
        end
    end
    diff_test_p(:,:,k)=pos_part+y_mean;
    diff_test_n(:,:,k)=neg_part+y_mean;
end
%计算检验标准
diff_test(:,:,:,1)=diff_test_p;
diff_test(:,:,:,2)=diff_test_n;
diff_test(:,:,:,3)=diff_test(:,:,:,1)-diff_test(:,:,:,2);
sigcrit=nan(sizy(1),sizy(2),3);
for k=1:3
    for k1=1:sizy(1)
        for k2=1:sizy(2)
         diff_test_sig=sort(diff_test(k1,k2,:,k),3);
         sigcrit1=diff_test_sig(1,1,floor(N*alpha0)+1);
         sigcrit2=diff_test_sig(1,1,floor(N*(1-alpha0))-1);
         if difference(k1,k2,k)<sigcrit1
             sigcrit(k1,k2,k)=sigcrit1;
         elseif difference(k1,k2,k)>sigcrit2
             sigcrit(k1,k2,k)=sigcrit2;
        end
    end
end

end
