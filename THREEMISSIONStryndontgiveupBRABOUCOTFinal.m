%% JAI GURU BABA
%% Hey Hori Hey Ram, Krishna Krishna, Bhagawanor Naam loi
% CHOOSING B MATRICES WITH OPTIMUM SOLUTION
% STEP 1 IS TO BUILD B MATRICES WITH NEGATIVE EIGEN VALUES
% STEP 2 IS TO FIND B MATRICES WITH OPTIMUM SOLUTION

% STEP 1 IS TO BUILD B MATRICES WITH NEGATIVE EIGEN VALUES
%Operation Completion Time Matrix
%*************************************************************************
%*************************************************************************
% ST AND CT FOR CALCULATING MATRIX A
tic
ST=[0];
M123=[5 5 5 5 5 -inf -inf -inf; 
      5 5 5 5 5    5 -inf -inf;
      5 5 5 5 5 -inf -inf -inf];
M321 = [0 0 0; 0 0 0; 0 0 0];  
[leni,brea]=size(M123);
naap=length(ST)+leni*brea;
M123_i=M123;
C_1=zeros(naap,naap);
C_1(:,:)=-inf;
for i=1 :naap
    C_1(i,i)=0;
    C_1(i,1)=0;    
end

[rows,cols]=size(M123);
for i=1:rows-2
      for j=1:2
          if i<2
        C_1(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_1(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

A_1=zeros(naap,naap);
A_1(:,:)=-inf;
for i=1 :naap
    A_1(i,i)=0;
    A_1(i,1)=0;    
end
for i=1:rows-2
      for j=3:cols
          if i<2
        A_1(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_1(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

MaxLAG1=-inf(naap,naap);

MaxLAG1(1,3)=-35;
MaxLAG1(1,4)=-35;
MaxLAG1(1,5)=-35;
MaxLAG1(1,6)=-35;
MaxLAG1(1,7)=-35;
%  LAG(1,6)=-5*380+60;LAG(1,11)=-5 *181+54;LAG(1,16)=-5 *239+55;
%  LAG(1,21)=-5 *195+37;
% if isempty(ST)~=1
% [row,col]=size(ST);
%     for j=1:col
%         A(j+1,1)=ST(1,j) ; 
%         
%     end   
%     
% end

C_2=zeros(naap,naap);
C_2(:,:)=-inf;
for i=1 :naap
    C_2(i,i)=0;
    C_2(i,1)=0;    
end

[rows,cols]=size(M123);
for i=2:rows-1
      for j=1:3
          if i<2
        C_2(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_2(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

A_2=zeros(naap,naap);
A_2(:,:)=-inf;
for i=1 :naap
    A_2(i,i)=0;
    A_2(i,1)=0;    
end
for i=2:rows-1
      for j=4:cols
          if i<2
        A_2(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_2(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

MaxLAG2=-inf(naap,naap);
MaxLAG2(1,11)=-40;
MaxLAG2(1,12)=-40;
MaxLAG2(1,13)=-40;
MaxLAG2(1,14)=-40;
MaxLAG2(1,15)=-40;
MaxLAG2(1,16)=-40;



C_3=zeros(naap,naap);
C_3(:,:)=-inf;
for i=1 :naap
    C_3(i,i)=0;
    C_3(i,1)=0;    
end

[rows,cols]=size(M123);
for i=3:rows
      for j=1:2
          if i<2
        C_3(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_3(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

A_3=zeros(naap,naap);
A_3(:,:)=-inf;
for i=1 :naap
    A_3(i,i)=0;
    A_3(i,1)=0;    
end
for i=3:rows
      for j=3:8
          if i<2
        A_3(j+cols*(i-1)+1+length(ST),j+1)=M123(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_3(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M123(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

MaxLAG3=-inf(naap,naap);
MaxLAG3(1,19)=-37;
MaxLAG3(1,20)=-37;
MaxLAG3(1,21)=-37;
MaxLAG3(1,22)=-37;
MaxLAG3(1,23)=-37;

A1=addition(C_1,A_1);
A_1=addition(A1,MaxLAG1); 

A2=addition(C_2,A_2);
A_2=addition(A2,MaxLAG2); 

A3=addition(C_3,A_3);
A_3=addition(A3,MaxLAG3); 

A4=addition(A_1,A_2);

%% So the final A matrix is

A=addition(A4,A_3);

% [cyvec, eigv, circ]=HowardA(A); 
% 
% eig=max(cyvec);

%% Calculating B matrix

T1=[ 5 5 5;    5 5 5;    5 5 5];
Iindex=[1 12 20; 4 10 18; 5 11 19];

[bapon, phukan] = size(T1); 

[sunu munu]= size(Iindex); 

for i=1:bapon 
    Ttemp(:,:,i)=perms(T1(i,:));
V2(:,:,i)=perms(Iindex(i,:));
end
   

%%Calculation for rp1, i.e. for segment A, we have 6 combinations, 

[nak,sah]=size(Iindex);
mac=nak;
N=sah;len=0;
i=1;ii=i;p=mac; storV=zeros(ii,N,mac);storT=zeros(ii,N,mac);
ne=0;
[nak,sah]=size(Iindex);
storgoodV1=[];storgoodT1=[];storelen=[]; 
for j=1:nak 
joBI=Iindex(j,:);

T=T1(j,:);ge=0;
[ goodV1,goodT1 ] = seeds( N,joBI,T,A); 
[rot,len]=size( goodV1);
storelen=[storelen;rot]; 
storgoodV1=[storgoodV1;goodV1];
storgoodT1=[storgoodT1 ;goodT1];
end
Ball=0;
for iml=1:length(storelen)
 Ball=Ball+ storelen(iml);
end
[bapa, bhanti, bablu]=size(V2); 
hor=25;ver=25;

W=zeros(hor,ver,bapa,bhanti);
K=-inf(hor,ver,bapa,bhanti);
W(:,:,:,:)=-inf;




% count=0;    
% for l=1:bapon    
%     dTtemp=Ttemp(:,:,l);
%     V=V2(:,:,l) ;
%    for i=1:bapa    % permutation of 1 op set
%        for j=2:phukan  % columns of T
%          for   k=1:j-1        
%             count=count+1  ;      
%     W(V(i,j)+1,V(i,k)+1,i,l)= dTtemp(i,k);   
%     
%          end        
%       end
%     end
% end

% [Bmax] = addition( W(:,:,1,3),W(:,:,1,3) );
% [BmaxF]=addition(Bmax,W(:,:,1,3));
% E=zeros(10,1);  
% A + B 
% APB=addition(A,W(:,:,4,3)); 
% [cyvec, eigv, circ]=HowardA(APB);   

SG1=storgoodV1;
ST1=storgoodT1;
NL1=storelen;

bre=1;
[koka aita]=size(SG1);
[amy ankit]=size(ST1); 
[riku siku]=size(A);
W1=-inf(riku, siku, koka, bre);
for l=1:bre  
    dTtemp=ST1(:,:,l);
    V=SG1(:,:,l) ;
   for n=1:koka    % permutation of 1 op set
       for j=2:ankit  % columns of T
         for   k=1:j-1        
                 
    W1(V(n,j)+1,V(n,k)+1,n,l)= dTtemp(n,k);    
         end        
      end
    end
end

 operations=leni*brea;  

 OW=[];NLc=0;

OW=W1; 
 NLc=NL1;
 ONLc=NL1;

for i=2:length(ONLc)
ONLc(i)=ONLc(i)+ONLc(i-1);
end

for i=2
NLc(i)=NLc(i)+NLc(i-1);
end

stseq=[];stW=[];Orggv1=[];
gooseq=[];

 Orggv1=storgoodV1;

 [hamsa hamsi]=size(storgoodV1); 

 counts=1;ceig=[];
 gooseq=[];gooW=[];storseq=[];

for c=1:length(NL1)-1
       count=0;
       toc
       
    for a=1:NLc(c)
    for b=ONLc(c)+1:ONLc(c+1)
  %###################################### 
        if counts==0  
   
    W1=gooW;
    storgoodV1=gooseq;
   for zz=2:length(NLc)
   NLc(zz)= NLc(zz)+NLc(zz-1) ;  
    end          
            gooseq=[];
%         
        end
    end   
% %        
% %        %########################################      
      counts=counts+1;
        sW=0;
        sW=addition(W1(:,:,a),OW(:,:,b));
        APB=addition( A,sW );
[cyvec, eigv, circ]=HowardA(APB);
eig=max(cyvec);
ceig=[ceig eig];
if eig==0
count=count+1;  
 st=[storgoodV1(a,:) Orggv1(b,:)];
    gooseq(count,:)=st; 
    gooW(:,:,count)=sW;
    st=[];
 end
  end
  end
   [nl,nw]=size(gooseq);
    NLc(c)=0;
    NLc(c+1)=nl;
    counts=0;
    [nib,mir]=size(gooseq);
    if mir==operations
        stseq=gooseq;
       stW=gooW;
    end

%% Shared resources among two AGVs

Iindexi=[2 21; 3 17];
T2=[5 5;
    5 5];

[damP ratan]=size(T2);

[kiku rini]=size(Iindexi); 

for i=1:damP 
    Ttemp1(:,:,i)=perms(T2(i,:)); 
    for p=1:kiku
        V3(:,:,p)=perms(Iindexi(p,:));
    end
end

[bhaiti,akon]=size(Iindexi);
mac=bhaiti;
N=akon;len=0;
i=1;ii=i;p=mac; storV1=zeros(ii,N,mac);storT1=zeros(ii,N,mac);
ne=0;
[bhaiti,akon]=size(Iindexi);
storgoodV2=[];storgoodT2=[];storelen1=[]; 
for j=1:bhaiti 
joBI=Iindexi(j,:);

Tm=T2(j,:);ge=0;
[ goodV1,goodT1 ] = seeds( N,joBI,T2,A); 
[rot,len]=size( goodV1);
storelen1=[storelen1;rot]; 
storgoodV2=[storgoodV2;goodV1];
storgoodT2=[storgoodT2 ;goodT1];
end
Ball1=0;
for iml=1:length(storelen1)
 Ball1=Ball1+ storelen1(iml);
end

SG2=storgoodV2;
ST2=storgoodT2;

bre1=1;
[vik prem]=size(SG2);
[kak ding]=size(ST2); 
[ravi kisan]=size(A);
NL2=storelen1;

W2=-inf(ravi, kisan, vik, bre1);
for f=1:bre1  
    dTtemp=ST2(:,:,f);
    V=SG2(:,:,f) ;
   for n=1:vik    % permutation of 1 op set
       for j=2:ding  % columns of T
         for   k=1:j-1        
                 
    W2(V(n,j)+1,V(n,k)+1,n,f)= dTtemp(n,k);    
         end        
      end
    end
end

operations1=leni*brea;  

 OW1=[];NLc1=0;

OW1=W2; 
 NLc1=NL2;
 ONLc1=NL2;

% for i=2:length(ONLc1)
% ONLc1(i)=ONLc1(i)+ONLc1(i-1);
% end
% 
% for i=2
% NLc1(i)=NLc1(i)+NLc1(i-1);
% end

stseq1=[];stW1=[];Orggv2=[];
gooseq1=[];

 Orggv2=storgoodV2;

 [loni moni]=size(storgoodV2); 

 counts1=1;ceig1=[];
 gooseq1=[];gooW1=[];storseq1=[];

for c=1:length(NL2)
       count1=0;
       toc
       
    for a=1:NLc1(c)
    for b=ONLc1(c):ONLc1(c)
  %###################################### 
        if counts1==0  
   
    W2=gooW1;
    storgoodV2=gooseq1;
   for zz=2:length(NLc1)
   NLc1(zz)= NLc1(zz)+NLc1(zz-1) ;  
    end          
            gooseq1=[];
%         
        end
    end   
% %        
% %        %########################################      
      counts1=counts1+1;
        sW1=0;
        sW1=addition(W2(:,:,a),OW1(:,:,b));
        APB1=addition( A,sW1 );
[cyvec, eigv, circ]=HowardA(APB1);
eig=max(cyvec);
ceig1=[ceig1 eig];
if eig==0
count1=count1+1;  
 st1=[storgoodV2(a,:) Orggv2(b,:)];
    gooseq1(count1,:)=st1; 
    gooW1(:,:,count1)=sW1;
    st1=[];
 end
  end
  end
   [nl1,nw1]=size(gooseq1);
    NLc1(c)=0;
    NLc1(c+1)=nl1;
    counts1=0;
    [nib1,mir1]=size(gooseq1);
    if mir1==operations1
        stseq1=gooseq1;
       stW1=gooW1;
    end

   
      [sdo,jto,de]=size(gooW);
      [tin,tim,tq]=size(gooW1);

  E=zeros(operations+1,1);Ti=0;T=0;

for soo=1:leni
      Ti=[Ti M123(soo,:)]; 
end

T3=Ti';
  CT=[];DT=[];Tf=[];    

for y=1:de
    for z=1:tq
  APB1=addition(A,gooW(:,:,y));
  APB3=addition(A,gooW1(:,:,z));
  APB4=addition(APB1,APB3);
 Akl  = powers( APB4,operations );
 D=multiplication( Akl,E );
 CT(:,y)=D;
 DT(:,y)=D+T3; 
    end
end
for s=1:de
    Tf(s)=max(DT(:,s));
end
[indx1 mintf]=min(Tf);

