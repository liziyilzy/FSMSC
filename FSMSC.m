load('PaviaU.mat');
I=paviaU;
I=double(I);
[A,B,C]=size(I);
% for i=1:C
%     I(:,:,i)=Normalize(I1(:,:,i));
% end
Z=reshape(I,A*B,C);
c=A*B;
K =9;
LC=1000;
Lm=1.3;
LM=2/(Lm-1);
iter=0;
epsm=1.0e-6;
[LN,LS]=size(Z);
LU(LC,LN)=0;
Dist(LC,LN)=0;
LP(LC,LS)=0;
LU0=rand(LC,LN);
LA=ones(LC,1);
LA=LA*sum(LU0);
LU0=LU0./LA;
while iter<1000
    iter=iter+1;
    %计算或更新聚类中心
    LUm=LU0.^Lm;
    LUm1=LUm*Z;
    LP=LUm1./(ones(LS,1)*sum(LUm'))';
    %更新划分矩阵U
    for i=1:LC
        for j=1:LN
            Dist(i,j)=fuzzydist(LP(i,:),Z(j,:));
        end
    end
    LU=1./(Dist.^Lm.*(ones(LC,1)*sum(Dist.^(-Lm))));
    if norm(LU-LU0,Inf)<epsm
        break
    end
    LU0=LU;
end
W=LU';
Scolumn =sum(W);
sco=diag(Scolumn);
sco=sco^(-1/2);
for ia=1:1000
    for ja=1:1000
        if sco(ia,ja)==inf;
            sco(ia,ja)=0;
        end
    end
end
W=W*sco;
[U,S,V]=svds(W,K);
for i1=1:K
    U1(:,i1)=U(:,i1);
end
Sumrow=sum(U1.*U1,2);
Sumrow=Sumrow.^0.5;
for ii=1:c
    for jj=1:K
        U2(ii,jj)=U1(ii,jj)/Sumrow(ii,1);
    end
end
[idx,ctrs] = kmeans(U2,K);%ctrs为中心点，idx为类别数(列矩阵)。
c1 = find(idx==1);
c2 = find(idx==2);
c3 = find(idx==3);
c4 = find(idx==4);
c5 = find(idx==5);
c6 = find(idx==6);
c7 = find(idx==7);
c8 = find(idx==8);
c9 = find(idx==9);
% c10 = find(idx==10);
% c11 = find(idx==11);
% c12 = find(idx==12);
% c13 = find(idx==13);
% c14 = find(idx==14);
% c15 = find(idx==15);
% c16 = find(idx==16);
Z1=zeros(c,3);
c1=c1';
[g,v]=size(c1);
for o=1:v
    Z1(c1(o),1)= 255;
    Z1(c1(o),2)= 222;
    Z1(c1(o),3)= 173;
end
c2=c2';
[u,t]=size(c2);
for s=1:t
    Z1(c2(s),1)=230;
    Z1(c2(s),2)=230;
    Z1(c2(s),3)=250;
end
c3=c3';
[u1,t1]=size(c3);
for s1=1:t1
    Z1(c3(s1),1)=105;
    Z1(c3(s1),2)=105;
    Z1(c3(s1),3)=105;
end
c4=c4';
[u2,t2]=size(c4);
for s2=1:t2
    Z(c4(s2),1)=0;
    Z(c4(s2),2)=0;
    Z(c4(s2),3)=128;
end
c5=c5';
[u3,t3]=size(c5);
for s3=1:t3
    Z(c5(s3),1)=135;
    Z(c5(s3),2)=206;
    Z(c5(s3),3)=235;
end
c6=c6';
[u4,t4]=size(c6);
for s4=1:t4
    Z(c6(s4),1)=0;
    Z(c6(s4),2)=255;
    Z(c6(s4),3)=255;
end
c7=c7';
[u5,t5]=size(c7);
for s5=1:t5
    Z(c7(s5),1)=0;
    Z(c7(s5),2)=100;
    Z(c7(s5),3)=0;
end
c8=c8';
[u6,t6]=size(c8);
for s6=1:t6
    Z(c8(s6),1)=0;
    Z(c8(s6),2)=255;
    Z(c8(s6),3)=127;
end
c9=c9';
[u7,t7]=size(c9);
for s7=1:t7
    Z(c9(s7),1)=154;
    Z(c9(s7),2)=205;
    Z(c9(s7),3)=50;
end
% c10=c10';
% [u8,t8]=size(c10);
% for s8=1:t8
%     Z(c10(s8),1)=255;
%     Z(c10(s8),2)=215;
%     Z(c10(s8),3)=0;
% end
% c11=c11';
% [u9,t9]=size(c11);
% for s9=1:t9
%     Z(c11(s9),1)=205;
%     Z(c11(s9),2)=92;
%     Z(c11(s9),3)=92;
% end
% c12=c12';
% [u10,t10]=size(c12);
% for s10=1:t10
%     Z(c12(s10),1)=139;
%     Z(c12(s10),2)=69;
%     Z(c12(s10),3)=19;
% end
% c13=c13';
% [u11,t11]=size(c13);
% for s11=1:t11
%     Z(c13(s11),1)=255;
%     Z(c13(s11),2)=140;
%     Z(c13(s11),3)=0;
% end
% c14=c14';
% [u12,t12]=size(c14);
% for s12=1:t12
%     Z(c14(s12),1)=255;
%     Z(c14(s12),2)=20;
%     Z(c14(s12),3)=147;
% end
% c15=c15';
% [u13,t13]=size(c15);
% for s13=1:t13
%     Z(c15(s13),1)=46;
%     Z(c15(s13),2)=139;
%     Z(c15(s13),3)=87;
% end
% c16=c16';
% [u14,t14]=size(c16);
% for s14=1:t14
%     Z(c16(s14),1)=0;
%     Z(c16(s14),2)=139;
%     Z(c16(s14),3)=139;
Y=reshape(Z1,A,B,3);
Y1=Y(:,:,1:3);
figure
imshow(Y1/255);
Y2=Y1/255;