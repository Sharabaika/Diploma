aMZ=1;              %0 - to calculate the array of neigbours, 1 - to read from disk
QW=0.7;             %relaxation parameters
QPsi=1.2;
QT=1.2;
NumbNodes=length(p);        %Number of nodes
NumbElements=length(t);     %Number of elements
Ra=1e8;                   %Rayleigh number
R_outer=1;                  %Outer radius
R_inner=0.95;               %Inner radius
T_outer=0;
T_inner=1;
Pr=10;
is_wall=zeros(1,NumbNodes);
Psi=zeros(1,NumbNodes);
T=zeros(1,NumbNodes);
W=zeros(1,NumbNodes); 
for i=1:NumbNodes
        %Determining whether a node belongs to the boundaries and the internal area
        if abs(sqrt(p(i,1)^2+p(i,2)^2)-R_inner) < 1e-10
            Psi(i)=0;
            T(i)=T_inner;
            is_wall(i)=1;       %Inner boundary of the annulus
        elseif abs(sqrt(p(i,1)^2+p(i,2)^2)-R_outer) < 1e-10
            Psi(i)=0;
            T(i)=T_outer;
            is_wall(i)=2;       %Outer boundary of the annulus
        else
            Psi(i)=0.0001;
            is_wall(i)=0;      %Nodes in the annulus
        end
end
%Construction of the array of neigbours: first column - number of neibour elements
%another columns - numbers of neigbour elements
if aMZ==0
        MZn=zeros(NumbNodes,9);
        for i=1:length(p)    %i - номер узла
            k=2;
            for j=1:length(t)    %j - номер элемента
                if i==t(j,1)||i==t(j,2)||i==t(j,3)
                    MZn(i,k)=j;
                    k=k+1;
                end
            end
            MZn(i,1)=k-2;
        end
        fmz=fopen('MZn_file.txt','wt');
        fwrite(fmz,MZn,'integer*4');
        fclose(fmz);
else        %Read from disk
        fmz=fopen('MZn_file.txt','rt' );
        MZn=fread(fmz,[NumbNodes 9],'integer*4');
        fclose(fmz);
end
S='Step MZn';
disp(S);
%End of construction of the array of neigbours

inaccPsi=zeros(1,NumbNodes);   %Zeroing of inaccuracies
inaccW=zeros(1,NumbNodes);
inaccT=zeros(1,NumbNodes);

st=1;   %Initial inaccuracy
n=0;    %Iteration counter
while st>1e-7     %Iteration loop
  sPsi=0; dPsi=0; sW=0; dW=0; sT=0; dT=0;
  %Node loop
  for i=1:NumbNodes
    aPsi0=0; aPsinb=0; aW0=0; aWnb=0; aT0=0; aTnb=0;
    S=0; I=0; SdTdx=0;
    PsiOld=Psi(i); WOld=W(i); TOld=T(i);
    %Neightbours loop
    for j=1:MZn(i,1)    
          nel=MZn(i,j+1);  %Number of the element
          n0=t(nel,1);n1=t(nel,2);n2=t(nel,3);  %Numbers of nodes
          if i==n1
              n1=n2;n2=n0;n0=i;
          end
          if i==n2
             n2=n1;n1=n0;n0=i;
          end                       %End of nodes renumeration
          x0=p(n0,1);y0=p(n0,2);
          x1=p(n1,1);y1=p(n1,2);
          x2=p(n2,1);y2=p(n2,2);
          x10=x1-x0;y01=y0-y1;
          x21=x2-x1;y12=y1-y2;
          x02=x0-x2;y20=y2-y0;
          x01=-x10;y10=-y01;
          x12=-x21;y21=-y12;
          x20=-x02;y02=-y20;
          Delta=x10*y20-x02*y01;
%%%%%%%%%%%%%%%%%%%%  Boundaries %%%%%%%%%%%%%%%%%%%%%%%%%
       if (is_wall(n0)==1) || (is_wall(n0)==2)
              if is_wall(n0) == 1
                  T(i)=1;
                  sina=-p(n0,1)/R_inner;
                  cosa=p(n0,2)/R_inner;
              elseif is_wall(n0)==2
                  T(i)=0;
                  sina=p(n0,1)/R_outer;
                  cosa=-p(n0,2)/R_outer;
              end
              X1=(x1-x0)*cosa+(y1-y0)*sina;
              Y1=(y1-y0)*cosa-(x1-x0)*sina;
              X2=(x2-x0)*cosa+(y2-y0)*sina;
              Y2=(y2-y0)*cosa-(x2-x0)*sina;
              if Y1==0
                  a=0;  b=2*(Psi(n2)-Psi(n0))/Y2^2;
              elseif Y2==0
                  a=0;  b=2*(Psi(n1)-Psi(n0))/Y1^2;
              else
                  DeltaGr=0.5*Y1*Y2*(X1*Y2-X2*Y1);
                  DeltaAGr=0.5*((Psi(n1)-Psi(n0))*Y2^2-(Psi(n2)-Psi(n0))*Y1^2);
                  DeltaBGr=(Psi(n2)-Psi(n0))*X1*Y1-(Psi(n1)-Psi(n0))*X2*Y2;
                  a=DeltaAGr/DeltaGr;
                  b=DeltaBGr/DeltaGr;
              end
              Xa=X1/2; Ya=Y1/2; Xb=X2/2; Yb=Y2/2;
              Xc=(X1+X2)/3; Yc=(Y1+Y2)/3;
              I=I+0.5*a*((Yb^2-Ya^2)-(Xb^2-Xa^2));
              I=I-0.5*b*((Ya+Yc)*(Xc-Xa)+(Yb+Yc)*(Xb-Xc));
              aWnb=aWnb+(7*W(n1)+7*W(n2))/216*Delta;
              aW0=aW0+11*Delta/108;
 %%%%%%%%%%%%%%%%%%%%  End dor boundaries %%%%%%%%%%%%%%%%%%%%%%%%%
       else         %Inner nodes
%%%%%%%%%%%%%%%%%%%%  Streamfunction %%%%%%%%%%%%%%%%%%%%%%%%%
          aPsi0=aPsi0+(x21^2+y12^2)/Delta;
          aPsinb=aPsinb+(Psi(n1)*(y20*y12+x02*x21)+Psi(n2)*(y01*y12+x10*x21))/Delta;
          S=S+(22*W(n0)+7*W(n1)+7*W(n2))*Delta/216;
%%%%%%%%%%%%%%%%%%%%  End for streamfunction %%%%%%%%%%%%%%%%%%%%%%%%%
          Apsi=Psi(n0)*y12+Psi(n1)*y20+Psi(n2)*y01;
          Bpsi=Psi(n0)*x21+Psi(n1)*x02+Psi(n2)*x10;
          Vx=Bpsi/Delta;
          Vy=-Apsi/Delta;
          Vel=sqrt(Vx^2+Vy^2);
          if abs(Vel)<1e-8
              cosa=1;
              sina=0;
          else
              cosa=Vx/Vel;
              sina=Vy/Vel;
          end
          xc=(x0+x1+x2)/3;
          yc=(y0+y1+y2)/3;  
          X0=(x0-xc)*cosa+(y0-yc)*sina;
          Y0=(y0-yc)*cosa-(x0-xc)*sina;
          X1=(x1-xc)*cosa+(y1-yc)*sina;
          Y1=(y1-yc)*cosa-(x1-xc)*sina;
          X2=(x2-xc)*cosa+(y2-yc)*sina;
          Y2=(y2-yc)*cosa-(x2-xc)*sina;
          X3=[X0,X1,X2];
          Xmax=max(X3); Xmin=min(X3);
          X12=X1-X2;Y12=Y1-Y2;
          Y20=Y2-Y0;Y01=Y0-Y1;  
          X10=X1-X0;X02=X0-X2;
%%%%%%%%%%%%%%%%%%%%  Temperature %%%%%%%%%%%%%%%%%%%%%%%%%
          abT=Vel*Y12*Y0/8+X12/2;
          acT=Vel*Y12/2;
          if abs(Vel)<1e-14
              k0T=0;
              k1T=0;
              k2T=0;
          elseif abs(Vel*(Xmax-Xmin))<1e-8
              DT=Vel*(X0*Y12+X1*Y20+X2*Y01);
              k0T=(abT*Vel*(X2-X1)+acT*(-Y12+Vel*((X1-Xmax)*Y2-(X2-Xmax)*Y1)))/DT;
              k1T=(abT*Vel*(X0-X2)+acT*(-Y20+Vel*((X2-Xmax)*Y0-(X0-Xmax)*Y2)))/DT;
              k2T=(abW*Vel*(X1-X0)+acT*(-Y01+Vel*((X0-Xmax)*Y1-(X1-Xmax)*Y0)))/DT;
          else
              E0T=exp(Vel*(X0-Xmax));
              E1T=exp(Vel*(X1-Xmax));
              E2T=exp(Vel*(X2-Xmax));
              DT=(E0T*Y12+E1T*Y20+E2T*Y01);
              k0T=(abT*(E2T-E1T)+acT*(E1T*Y2-E2T*Y1))/DT;
              k1T=(abT*(E0T-E2T)+acT*(E2T*Y0-E0T*Y2))/DT;
              k2T=(abT*(E1T-E0T)+acT*(E0T*Y1-E1T*Y0))/DT;
          end
          AT=T(n0)*y12+T(n1)*y20+T(n2)*y01;
          SdTdx=SdTdx+AT/6;
          aT0=aT0+k0T;
          aTnb=aTnb+k1T*T(n1)+k2T*T(n2);
%%%%%%%%%%%%%%%%%%%%  End for temperature %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Vorticity for inner nodes %%%%%%%%%%%%%%%%%%%%%%%%%
              VelPr=Vel/Pr;
              abW=VelPr*Y12*Y0/8+X12/2;
              acW=VelPr*Y12/2;
              if abs(VelPr)<1e-15
                  k0W=0;
                  k1W=0;
                  k2W=0;
              elseif abs(VelPr*(Xmax-Xmin))<1e-8
                    DW=VelPr*(X0*Y12+X1*Y20+X2*Y01);
                    k0W=(abW*VelPr*(X2-X1)+acW*(-Y12+VelPr*((X1-Xmax)*Y2-(X2-Xmax)*Y1)))/DW;
                    k1W=(abW*VelPr*(X0-X2)+acW*(-Y20+VelPr*((X2-Xmax)*Y0-(X0-Xmax)*Y2)))/DW;
                    k2W=(abW*VelPr*(X1-X0)+acW*(-Y01+VelPr*((X0-Xmax)*Y1-(X1-Xmax)*Y0)))/DW;
              else
                    E0=exp(VelPr*(X0-Xmax));
                    E1=exp(VelPr*(X1-Xmax));
                    E2=exp(VelPr*(X2-Xmax));
                    DW=(E0*Y12+E1*Y20+E2*Y01);
                    k0W=(abW*(E2-E1)+acW*(E1*Y2-E2*Y1))/DW;
                    k1W=(abW*(E0-E2)+acW*(E2*Y0-E0*Y2))/DW;
                    k2W=(abW*(E1-E0)+acW*(E0*Y1-E1*Y0))/DW;
              end
              aW0=aW0+k0W;
              aWnb=aWnb+k1W*W(n1)+k2W*W(n2);
       end
    end %End of neighbours loop
        if is_wall(i)==0
            Psi(i)=(-aPsinb+S)/aPsi0;
            Psi(i)=Psi(i)*QPsi+PsiOld*(1-QPsi);
            inaccPsi(i)=abs(PsiOld-Psi(i));
            W(i)=-(aWnb+Ra*SdTdx)/aW0;
            W(i)=W(i)*QW+WOld*(1-QW);
            inaccW(i)=abs(WOld-W(i));
            T(i)=-aTnb/aT0;
            T(i)=T(i)*QT+TOld*(1-QT);
            inaccT(i)=abs(TOld-T(i));
        else
            Psi(i)=0;
            W(i)=-(I+aWnb)/aW0;
            W(i)=W(i)*QW+WOld*(1-QW);
            inaccW(i)=abs(WOld-W(i));
            if is_wall(i)==1
                T(i)=T_inner;
            elseif is_wall(i)==2
                T(i)=T_outer;
            end
        end
    end       %End of of nodes loop
    if sum(sum(abs(Psi)))~=0||sum(sum(abs(W)))~=0
        stp=sum(sum(inaccPsi))/sum(sum(abs(Psi)))/QPsi;
        stw=sum(sum(inaccW))/sum(sum(abs(W)))/QW;
        stt=sum(sum(inaccT))/sum(sum(abs(T)))/QT;
        st=max(stp,stw);
        st=max(st,stt);
    else
        st=1;
    end
    n=n+1;
    if rem(n,100)==0
      disp(n);disp(st);
    end
end   %End of iteration loop
%picture
figure(3);
axis equal
Tmin=min(Psi);
DeltaT=max(Psi)-Tmin; 
C=zeros(1,3);
for i=1:length(t)
    n0=t(i,1); n1=t(i,2); n2=t(i,3);
    XX=[p(n0,1) p(n1,1) p(n2,1)];
    YY=[p(n0,2) p(n1,2) p(n2,2)];
    TT=((Psi(n0)+Psi(n1)+Psi(n2))/3-Tmin)/DeltaT;
    if TT<1/3
        C(1)=0; C(2)=3*TT; C(3)=1;
    elseif (TT>=1/3)&&(TT<2/3)
        C(1)=3*(TT-1/3); C(2)=1; C(3)=1-3*(TT-1/3);
    elseif (TT>=2/3)&&(TT<=1)
        C(1)=1; C(2)=1-3*(TT-2/3); C(3)=0;
    end
    fill(XX,YY,C,'LineStyle','none');
    hold on;
end
figure(4);
Tmin=min(T);
DeltaT=max(T)-Tmin; 
C=zeros(1,3);
for i=1:length(t)
    n0=t(i,1); n1=t(i,2); n2=t(i,3);
    XX=[p(n0,1) p(n1,1) p(n2,1)];
    YY=[p(n0,2) p(n1,2) p(n2,2)];
    TT=((T(n0)+T(n1)+T(n2))/3-Tmin)/DeltaT;
    if TT<1/3
        C(1)=0; C(2)=3*TT; C(3)=1;
    elseif (TT>=1/3)&&(TT<2/3)
        C(1)=3*(TT-1/3); C(2)=1; C(3)=1-3*(TT-1/3);
    elseif (TT>=2/3)&&(TT<=1)
        C(1)=1; C(2)=1-3*(TT-2/3); C(3)=0;
    end
    fill(XX,YY,C,'LineStyle','none');
    hold on;
end
