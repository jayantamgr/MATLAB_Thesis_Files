%%Hey Hori Hey Ram, Bhagawanor Namm loi, Krishna Krishna
%%A new test bed of 13 segments, B and G segments are priviledged based
%%segments. It is priviledged to only AGV with full skills to reduce
%%congestion. Anyway it has to travel the wwhole test bed. But segments,
%%A I M K F D will be used only by the AGVs with skills having Vertical
%%transfer and 15 kg load. Similarly M L H E C J segments will be used by
%%the AGVs with skills of horiuontal tranfer and 15 kg load.

%%Conditions---
%%St1 requests for an AGV for a transfer to St5, an AGV which is AGV 2 is
%%parked in St3 travels to S1 for loading the materials and finals goes to
%%St5 for unloading. That makes a mission, Mission 1. Mission 2 is where an
%%AGV is requested by St2. A parked AGV from St5 goes to St2 for loading
%%and has to go to St4 for unloading. Mission 3 is where an AGV is
%%requested by St5, AGV5 parked at St2 goes to St1 for loading and has to
%%delivar to St3. Mission 4 is an AGV from St3 has the request, so an
%%already parked AGV in St3 loads the materials at St3 and goes to station
%%1 for unloading(Note- Two AGVs starts from same station at a time, a segment
%%can be used by one AGV at a time, so we have minimum lag constraints here).
%%Mission 5 is requested by St4, so a parked AGV in St4 goes to St5 for
%%loading and goes to St2 for unloading(same constraints apply here also like
%%Mission 4).

%%Defining the segments
%%AGV speed is 3m/s, so we have to calculate the time that an AGV will
%%take in all the segments.

A=7; %%Station1
B=12;
C=7; %%Station2
D=11;
E1=12;
F=7; %%Station3
G=11; 
H=7; %%Station4
I=8;
J=6;
K=8;
L=7;
M=7; %%Station5

tic
ST=[0];

M12=[D  A  I    M -inf -inf -inf -inf -inf -inf -inf -inf -inf;
%      J C E1    H -inf -inf -inf -inf -inf -inf -inf -inf -inf;
     E1 H  L    M -inf -inf -inf -inf -inf -inf -inf -inf -inf;
     B  I  M    K    F -inf -inf -inf -inf -inf -inf -inf -inf;
%      F D A -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf;
     C  E1 H     -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf;
     L  M  J    C -inf -inf -inf -inf -inf -inf -inf -inf -inf;
     ];            
    
[leni,brea]=size(M12);
naap=length(ST)+leni*brea;
Mi1=M12;
Mi2=M12;
C_xot=zeros(naap,naap);
C_xot(:,:)=-inf;

for i=1 :naap
   C_xot(i,i)=0;
    C_xot(i,1)=0;    
end

[rows,cols]=size(M12);

for i=1:rows-4
      for j=1:2
          if i<2
        C_xot(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_xot(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

 
A_xot=zeros(naap,naap);
A_xot(:,:)=-inf;
for i=1 :naap
   A_xot(i,i)=0;
    A_xot(i,1)=0;    
end
Mis1=M12;


for i=1:rows-4
      for j=3:cols
          if i<2
        A_xot(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_xot(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end 


MaxLAG=-inf(naap,naap);
% MaxLAG(1,25)=-5;

MaxLAG(1,3)=-90;
MaxLAG(1,4)=-90;
MaxLAG(1,5)=-90;
MaxLAG(1,6)=-90;

% [rr,cc]=size(MaxLAG);
% for i = 1:1
%     for j = 1+10:cc-8;
%       % MaxLAG(1,i)=-25;
%        MaxLAG(i,j)=-25;
%     end
%     
% end

A1=addition(C_xot, A_xot);
A2=addition(A1,MaxLAG);

C_xot1=zeros(naap,naap);
C_xot1(:,:)=-inf;

for i=1 :naap
   C_xot1(i,i)=0;
    C_xot1(i,1)=0;    
end

[rows,cols]=size(M12);
for i=2:rows-3
      for j=1:2
          if i<2
        C_xot1(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_xot1(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end
 
A_xot1=zeros(naap,naap);
A_xot1(:,:)=-inf;

for i=1 :naap
    A_xot1(i,i)=0;
    A_xot1(i,1)=0;    
end
[rows,cols]=size(M12);
for i=2:rows-3
      for j=3:cols
          if i<2
        A_xot1(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_xot1(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

      

 

MaxLAG1=-inf(naap,naap);
% MaxLAG1(1,25)=-30;

MaxLAG1(1,16)=-80;
MaxLAG1(1,17)=-80;
MaxLAG1(1,18)=-80;
MaxLAG1(1,19)=-80;





A3=addition(C_xot1,A_xot1);
A4=addition(A3, MaxLAG1);

A5=addition(A2,A4);


C_xot2=zeros(naap,naap);
C_xot2(:,:)=-inf;

for i=1 :naap
   C_xot2(i,i)=0;
    C_xot2(i,1)=0;    
end

[rows,cols]=size(M12);

for i=3:rows-2
      for j=1:2
          if i<2
        C_xot2(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_xot2(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

 
A_xot2=zeros(naap,naap);
A_xot2(:,:)=-inf;
for i=1 :naap
   A_xot2(i,i)=0;
    A_xot2(i,1)=0;    
end
 

for i=3:rows-2
      for j=3:cols
          if i<2
        A_xot2(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_xot2(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end
      

MaxLAG2=-inf(naap,naap);
% MaxLAG(1,25)=-5;

MaxLAG2(1,29)=-85;
MaxLAG2(1,30)=-85;
MaxLAG2(1,31)=-85;
MaxLAG2(1,32)=-85;
MaxLAG2(1,33)=-85;


A6=addition(C_xot2, A_xot2);
A7=addition(A6,MaxLAG2);
A8=addition(A5,A7);


C_xot3=zeros(naap,naap);
C_xot3(:,:)=-inf;

for i=1 :naap
   C_xot3(i,i)=0;
    C_xot3(i,1)=0;    
end

 
[rows,cols]=size(M12);
for i=4:rows-1
      for j=1:2
          if i<2
        C_xot3(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_xot3(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end
 
A_xot3=zeros(naap,naap);
A_xot3(:,:)=-inf;

for i=1 :naap
    A_xot3(i,i)=0;
    A_xot3(i,1)=0;    
end
[rows,cols]=size(M12);
for i=4:rows-1
      for j=3:cols
          if i<2
        A_xot3(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_xot3(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

      

 

MaxLAG3=-inf(naap,naap);
% MaxLAG1(1,25)=-30;

MaxLAG3(1,42)=-70;
MaxLAG3(1,43)=-70;
MaxLAG3(1,44)=-70;
%%Minimum time lag
MaxLAG3(42,41)=12;

%%Minimum time lag for this mission
%MaxLAG3(54,15)=15;

A9=addition(C_xot3,A_xot3);
A10=addition(A9, MaxLAG3);

A11=addition(A8,A10);


C_xot4=zeros(naap,naap);
C_xot4(:,:)=-inf;

for i=1 :naap
   C_xot4(i,i)=0;
    C_xot4(i,1)=0;    
end
  
[rows,cols]=size(M12);
for i=5:rows
      for j=1:3
          if i<2
        C_xot4(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        C_xot4(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end
 
A_xot4=zeros(naap,naap);
A_xot4(:,:)=-inf;

for i=1 :naap
    A_xot4(i,i)=0;
    A_xot4(i,1)=0;    
end
[rows,cols]=size(M12);
for i=5:rows
      for j=4:cols
          if i<2
        A_xot4(j+cols*(i-1)+1+length(ST),j+1)=M12(i,j) ; 
         %A(j+1,j+cols*(i-1)+1+length(ST))=-SCT(i,j) ;
          end
        if i>=2             
        A_xot4(j+cols*(i-1)+1+length(ST),j+cols*(i-1)+length(ST))=M12(i,j);  
        %A(j+cols*(i-1)+length(ST),j+cols*(i-1)+1+length(ST))=-SCT(i,j); 
        end         
      end 
end

      

 

MaxLAG4=-inf(naap,naap);
% MaxLAG1(1,25)=-30;

MaxLAG4(1,55)=-55;
MaxLAG4(1,56)=-55;
MaxLAG4(1,57)=-55;
MaxLAG4(1,58)=-55;



A12=addition(C_xot4,A_xot4);
A13=addition(A12, MaxLAG4);

A=addition(A11,A13);

[cyvec, eigv, circ]=HowardA(A); 

eig=max(cyvec);

%% Shared resources among two AGVs

Iindex3=[3 28; 14 41; 15 42; 16 53; 40 56];
T3=[8 8;
   12 12;
     7 7;
     7 7;
     7 7];

[damP1 ratan1]=size(T3);

[kiku1 rini1]=size(Iindex3); 

for i=1:damP1 
    Ttemp3(:,:,i)=perms(T3(i,:));
    for p=1:kiku1
        V4(:,:,p)=perms(Iindex3(p,:));
    end
end

[bhaiti1,akon1]=size(Iindex3);
mac=bhaiti1;
N=akon1;len=0;
i=1;ii=i;p=mac; storV3=zeros(ii,N,mac);storT3=zeros(ii,N,mac);
ne=0;
[bhaiti1,akon1]=size(Iindex3);
storgoodV3=[];storgoodT3=[];storelen3=[]; 
for j=1:bhaiti1 
joBI=Iindex3(j,:);

Tm=T3(j,:);ge=0;
[ goodV1,goodT1 ] = seeds( N,joBI,T3,A); 
[rot,len]=size( goodV1);
storelen3=[storelen3;rot]; 
storgoodV3=[storgoodV3;goodV1];
storgoodT3=[storgoodT3 ;goodT1];
end
Ball3=0;
for iml=1:length(storelen3)
 Ball3=Ball3+ storelen3(iml);
end

SG3=storgoodV3;
ST3=storgoodT3;

bre3=1;
[vik1 prem1]=size(SG3);
[kak1 ding1]=size(ST3); 
[ravi1 kisan1]=size(A);
NL3=storelen3;

W3=-inf(ravi1, kisan1, vik1, bre3);
for f=1:bre3  
    dTtemp3=ST3(:,:,f);
    V3=SG3(:,:,f) ;
   for n=1:vik1    % permutation of 1 op set
       for j=2:ding1  % columns of T
         for   k=1:j-1        
                 
    W3(V3(n,j)+1,V3(n,k)+1,n,f)= dTtemp3(n,k);    
         end        
      end
    end
end

operations3=leni*brea;  

 OW3=[];NLc3=0;

OW3=W3; 
 NLc3=NL3;
 ONLc3=NL3;

% for i=2:length(ONLc1)
% ONLc1(i)=ONLc1(i)+ONLc1(i-1);
% end
% 
% for i=2
% NLc1(i)=NLc1(i)+NLc1(i-1);
% end

stseq3=[];stW3=[];Orggv3=[];
gooseq3=[];

 Orggv3=storgoodV3;

 [loni1 moni1]=size(storgoodV3); 

 counts3=1;ceig3=[];
 gooseq3=[];gooW3=[];storseq3=[];

for c=1:length(NL3)
       count3=0;
       toc
       
    for a=1:NLc3(c)
    for b=1:ONLc3(c)
  %###################################### 
        if counts3==0  
   
    W3=gooW3;
    storgoodV3=gooseq3;
   for zz=2:length(NLc3)
   NLc3(zz)= NLc3(zz)+NLc3(zz-1) ;  
    end          
            gooseq3=[];
%         
        end
    end   
% %        
% %        %########################################      
      counts3=counts3+1;
        sW3=0;
        sW3=addition(W3(:,:,a),OW3(:,:,b));
        APB33=addition( A,sW3 );
[cyvec, eigv, circ]=HowardA(APB33);
eig=max(cyvec);
ceig3=[ceig3 eig];
if eig==0
count3=count3+1;  
 st3=[storgoodV3(a,:) Orggv3(b,:)];
    gooseq3(count3,:)=st3; 
    gooW3(:,:,count3)=sW3;
    st3=[];
 end
  end
  end
   [nl3,nw3]=size(gooseq3);
    NLc3(c)=0;
    NLc3(c+1)=nl3;
    counts3=0;
    [nib3,mir3]=size(gooseq3);
    if mir3==operations3
        stseq3=gooseq3;
       stW3=gooW3;
    end

%% B matrices for booking the AGV, which AGV will go 

%%T1 is the CT of the shared segments by the AGVs
T1=[7 7 7 7];

%%Iindex is the index matrix of the segments, which transition is used and when and by whom
Iindex=[4 17 29 54];

[bapon, phukan] = size(T1); 

[sunu munu]= size(Iindex); 

for i=1:bapon 
    Ttemp(:,:,i)=perms(T1(i,:));
V2(:,:,i)=perms(Iindex(i,:));
end
   
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

% for i=1:length(ONLc)
% ONLc(i)=ONLc(i)+ONLc(i);
% end
% 
% for i=1
% NLc(i)=NLc(i)+NLc(i);
% end

stseq=[];stW=[];Orggv1=[];
gooseq=[];

 Orggv1=storgoodV1;

 [hamsa hamsi]=size(storgoodV1); 

 counts=1;ceig=[];
 gooseq=[];gooW=[];storseq=[];

 for c=1:length(NL1)
     count=0;
        toc
       
    for a=1:NLc(c)
    for b=1:ONLc(c)
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
    
[sdo,jto,de]=size(gooW3);
      [tin,tim,tq]=size(gooW);

  E=zeros(operations+1,1);Ti=0;T=0;

for soo=1:leni
      Ti=[Ti M12(soo,:)]; 
end

T3=Ti';
  COT=[];DOT=[];Tf=[];    

for y=1:de
    for z=1:tq
  APB3=addition(A,gooW(:,:,z));
  APB2=addition(A,gooW3(:,:,y));
  APB4=addition(APB2,APB3);
 Akl  = powers( APB4,operations );
 D=multiplication( Akl,E );
 COT(:,y)=D;
 DOT(:,y)=D+T; 
    end
end
for s=1:de
    Tf(s)=max(DOT(:,s));
end
[indx1 mintf]=min(Tf);



