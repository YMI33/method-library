function [difference sigcrit]=wcomposite_para(x,y,alpha0)
%�����Ȩ�ϳɷ����ĺ���������˫�����ؿ������,���а汾
%������ʽ [WCA significance]=wcomposite(index,field,alpha0)
%����indexΪһάָ����fieldΪ��Ȩ�ϳɷ�������������������Ϊ��ά����������άΪʱ��ά��alpha0��ͳ�Ƽ����׼����0.95��ʾ95%ͳ�Ƽ���
%�����WCAΪ��ά������ǰ��ά�������field��ͬ������άΪ��ͬ�쳣��������쳣���ΪWCA(:,:,1),���쳣ΪWCA(:,:,2)���ϳɲ�ΪWCA(:,:,3)
%significanceΪ������������significance(:,:,1)Ϊ�����Ͻ磬significance(:,:,2)Ϊ�����½�,WCAֵ�����Ͻ������½�ļ�Ϊͨ������
sizlen=size(x);
if sizlen(1)==1
    x=x';
    sizlen=sizlen(2);
else
    sizlen=sizlen(1);
end

%˫�����
alpha0=alpha0/2;

sizy=size(y);
siznum=3;
%����ϳɽ��
x_anomal=x-mean(x);
y_mean=mean(y,siznum);
y_anomal=y-repmat(y_mean,[ones(1,siznum-1) sizy(siznum)]);
p=find(x_anomal>0);
n=find(x_anomal<0);

%difference����ά�ֱ�Ϊ��1.positive part; 2.negative part; 3. differnece
difference=nan(sizy(1),sizy(2),3);
for k1=1:sizy(1)
    for k2=1:sizy(2)
      difference(k1,k2,1)=sum(x_anomal(p).*squeeze(y_anomal(k1,k2,p)))/sum(abs(x_anomal(p)))+y_mean(k1,k2);
      difference(k1,k2,2)=-sum(x_anomal(n).*squeeze(y_anomal(k1,k2,n)))/sum(abs(x_anomal(n)))+y_mean(k1,k2);
    end
end
difference(:,:,3)=difference(:,:,1)-difference(:,:,2);

%����ͳ�Ƽ���
N=200;%���ؿ�������������
diff_test=nan(sizy(1),sizy(2),N,3);
%����������еļ�Ȩ�ϳɷ������
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
%��������׼
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