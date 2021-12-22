aMZ=0;              %0 - считать массив связи, 1 - читать с диска
Nx=80;             %Кол-во отрезков по азимуту
hx=1/Nx;
Ny=80;              %Кол-во отрезков по радиусу
hy=1/Ny;
QW=0.4;             %Пар-ры релаксации
QPsi=1.2;
NumbNodes=(Nx+1)*(Ny+1);    %Кол-во узлов
NumbElements=2*Nx*Ny;       %Кол-во элементов
u0=1;         %Скорость крышки
Re = 10000.0; %Число Рейнольдса
%вычисляем координаты узлов и пр.
for i=1:Nx+1
    x(i)=hx*(i-1);
end
for j=1:Ny+1
    y(j)=hy*(j-1);
end
p=zeros(2,NumbNodes);
is_wall=zeros(1,NumbNodes);
Psi=zeros(1,NumbNodes);
for i=1:Nx+1
    for j=1:Ny+1
        ii=(j-1)*(Nx+1)+i;          %Индекс одномерных массивов
        p(1,ii)=x(i); %Координата х
        p(2,ii)=y(j); %Координата у
        %Определение принадлежности узла границам и внутренней области
        if i==1
            Psi(ii)=0;
            if j==1
                is_wall(ii)=5;  %Левый нижний угол
            elseif j==Ny+1
                is_wall(ii)=7;  %Левый верхний угол
            else
                is_wall(ii)=1;     %Левая стенка
            end
        elseif j==1
            Psi(ii)=0;
            if i==Nx+1
                is_wall(ii)=6;  %Правый нижний угол
            else
                is_wall(ii)=2;      %Нижняя стенка
            end
        elseif j==Ny+1
            Psi(ii)=0;
            if i==Nx+1
                is_wall(ii)=8;     %Правый верхний угол
            else
                is_wall(ii)=3;      %Верхняя стенка
            end
        elseif i==Nx+1
            Psi(ii)=0;
            is_wall(ii)=4;       %Правая стенка
        else
            Psi(ii)=0.0001*sin(x(i)*pi)*sin(y(j)*pi);
            is_wall(ii)=0;      %Внутренние точки
        end
    end
end
%Конец определения координат, скоростей и пр.
%Построение треугольных элементов
t=zeros(NumbElements,3);
for j=1:Ny
    for i=1:Nx
            nel=2*(i-1)+1+(j-1)*Nx*2;
            n1=i+(j-1)*(Nx+1);
            n3=n1+Nx+1;
            n2=n3+1;
            t(nel,1)=n1;
            t(nel,2)=n2;
            t(nel,3)=n3;
            nel=nel+1;
            n2=n1+1;
            n3=n2+Nx+1;
            t(nel,1)=n1;
            t(nel,2)=n2;
            t(nel,3)=n3;
    end
end
%Конец построения треугольных элементов
%Картинка сетки
% figure(1)
% axis equal
% hold on
% for i=1:NumbElements
%     xxx(1)=p(1,t(i,1));
%     xxx(2)=p(1,t(i,2));
%     xxx(3)=p(1,t(i,3));
%     xxx(4)=p(1,t(i,1));
%     yyy(1)=p(2,t(i,1));
%     yyy(2)=p(2,t(i,2));
%     yyy(3)=p(2,t(i,3));
%     yyy(4)=p(2,t(i,1));
%     plot(xxx,yyy,'Color',[0 0 0])
% end
%Создание массива связи: первый столбец - кол-во элтов-соседей,
%прочие - номера соседних с i-тым узлом (i - номер строки) элементов
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
else
        fmz=fopen('MZn_file.txt','rt' );
        MZn=fread(fmz,[NumbNodes 9],'integer*4');
        fclose(fmz);
end
S='Step MZn';
disp(S);
%Конец постоения массива связи
%Вычисление температуры
inaccPsi=zeros(1,NumbNodes);
inaccW=zeros(1,NumbNodes);
W=zeros(1,NumbNodes); 
st=1; n=0;
while st>1e-7
  sPsi=0; dPsi=0; sW=0; dW=0;
  %Цикл по узлам
  for i=1:length(p) % i цикл по узлам 
    aPsi0=0; aPsinb=0; aW0=0; aWnb=0; S=0; I=0;
    PsiOld=Psi(i); WOld=W(i);
        %Цикл по соседним элементам
        for j=1:MZn(i,1)    % цикл по эл-там вокруг узла
          nel=MZn(i,j+1);  %номер соседнего эл-та
          n0=t(nel,1);n1=t(nel,2);n2=t(nel,3);
          if i==n1
              n1=n2;n2=n0;n0=i;
          end
          if i==n2
             n2=n1;n1=n0;n0=i;
          end
          x0=p(1,n0);y0=p(2,n0);
          x1=p(1,n1);y1=p(2,n1);
          x2=p(1,n2);y2=p(2,n2);
          x10=x1-x0;y01=y0-y1;
          x21=x2-x1;y12=y1-y2;
          x02=x0-x2;y20=y2-y0;
          x01=-x10;y10=-y01;
          x12=-x21;y21=-y12;
          x20=-x02;y02=-y20;
          Delta=x10*y20-x02*y01;
          if is_wall(n0)~=0
              if is_wall(n0)==1
                  sina=-1; cosa=0;
              elseif is_wall(n0)==2
                  sina=0; cosa=1;
              elseif is_wall(n0)==5
                  if is_wall(n2)==1
                      sina=-1; cosa=0;
                  else
                      sina=0; cosa=1;
                  end
              elseif is_wall(n0)==6
                  if is_wall(n1)==4
                      sina=1; cosa=0;
                  else
                      sina=0; cosa=1;
                  end
              elseif is_wall(n0)==3
                  sina=0; cosa=-1;
              elseif is_wall(n0)==7
                  if is_wall(n1)==1
                      sina=-1; cosa=0;
                  else
                      sina=0; cosa=-1;
                  end
              elseif is_wall(n0)==8
                  if is_wall(n2)==4
                      sina=1; cosa=0;
                  else
                      sina=0; cosa=-1;
                  end
              elseif is_wall(n0)==4
                  sina=1; cosa=0;
              end
              if is_wall(n0)==3 ||is_wall(n0)==7||is_wall(n0)==8
                  U=u0;
              else
                  U=0;
              end
              X1=(x1-x0)*cosa+(y1-y0)*sina;
              Y1=(y1-y0)*cosa-(x1-x0)*sina;
              X2=(x2-x0)*cosa+(y2-y0)*sina;
              Y2=(y2-y0)*cosa-(x2-x0)*sina;
              U=U*cosa;   %U изначально направлено вдоль оси х
              if is_wall(n0)==1||is_wall(n0)==2||is_wall(n0)==3||is_wall(n0)==4
                  if Y1==0
                      a=0;  b=2*(Psi(n2)-Psi(n0)-U*Y2)/Y2^2;
                  elseif Y2==0
                      a=0;  b=2*(Psi(n1)-Psi(n0)-U*Y1)/Y1^2;
                  else
                      DeltaGr=0.5*Y1*Y2*(X1*Y2-X2*Y1);
                      DeltaAGr=0.5*((Psi(n1)-Psi(n0)-U*Y1)*Y2^2-(Psi(n2)-Psi(n0)-U*Y2)*Y1^2);
                      DeltaBGr=(Psi(n2)-Psi(n0)-U*Y2)*X1*Y1-(Psi(n1)-Psi(n0)-U*Y1)*X2*Y2;
                      a=DeltaAGr/DeltaGr;
                      b=DeltaBGr/DeltaGr;
                  end
              elseif is_wall(n0)==7||is_wall(n0)==6
                  if is_wall(n1)==4||is_wall(n1)==1
                      a=0; b=2*(Psi(n2)-Psi(n0)-U*Y2)/Y2^2;
                  elseif is_wall(n2)==2||is_wall(n2)==3
                      a=0; b=2*(Psi(n1)-Psi(n0)-U*Y1)/Y1^2;
                  end
              elseif is_wall(n0)==5||is_wall(n0)==8
                  if is_wall(n1)==2||is_wall(n1)==3
                      a=0; b=2*(Psi(n2)-Psi(n0)-U*Y2)/Y2^2;
                  elseif is_wall(n2)==1||is_wall(n2)==4
                      a=0; b=2*(Psi(n1)-Psi(n0)-U*Y1)/Y1^2;
                  end
              end
              Xa=X1/2; Ya=Y1/2; Xb=X2/2; Yb=Y2/2;
              Xc=(X1+X2)/3; Yc=(Y1+Y2)/3;
              I=I-U*(Xb-Xa)+0.5*a*((Yb^2-Ya^2)-(Xb^2-Xa^2));
              I=I-0.5*b*((Ya+Yc)*(Xc-Xa)+(Yb+Yc)*(Xb-Xc));
              if is_wall(n1)==3||is_wall(n1)==7||is_wall(n1)==8
                  I=I+a*Ya^2/2-(U*Xa+a*Xa^2/2+b*Xa*Ya/2);
              elseif is_wall(n2)==3||is_wall(n2)==7||is_wall(n2)==8
                  I=I-a*Yb^2/2+(U*Xb+a*Xb^2/2+b*Xb*Yb/2);
              end
              aWnb=aWnb+(7*W(n1)+7*W(n2))/216*Delta;
              aW0=aW0+11*Delta/108;
          else
            Apsi=Psi(n0)*y12+Psi(n1)*y20+Psi(n2)*y01;
            Bpsi=Psi(n0)*x21+Psi(n1)*x02+Psi(n2)*x10;
            Vx=Bpsi/Delta;
            Vy=-Apsi/Delta;
            Vel=sqrt(Vx^2+Vy^2);
            cosa=Vx/Vel;
            sina=Vy/Vel;
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
            abW=Vel*Re*Y12*Y0/8+X12/2;
            acW=Vel*Re*Y12/2;
            ReVel=Re*Vel;
            if ReVel*(Xmax-Xmin)<1e-8
                DW=ReVel*(X0*Y12+X1*Y20+X2*Y01);
                k0W=(abW*ReVel*(X2-X1)+acW*(-Y12+ReVel*((X1-Xmax)*Y2-(X2-Xmax)*Y1)))/DW;
                k1W=(abW*ReVel*(X0-X2)+acW*(-Y20+ReVel*((X2-Xmax)*Y0-(X0-Xmax)*Y2)))/DW;
                k2W=(abW*ReVel*(X1-X0)+acW*(-Y01+ReVel*((X0-Xmax)*Y1-(X1-Xmax)*Y0)))/DW;
            else
%                 if ReVel*(Xmin-Xmax)<-10
%                     if abs(Xmax-X0)<1e-10
%                         DW=Y12;
%                         E0=1; E1=0; E2=0;
%                     elseif abs(Xmax-X1)<1e-10
%                         DW=Y20;
%                         E0=0; E1=1; E2=0;
%                     elseif abs(Xmax-X2)<1e-10
%                         DW=Y01;
%                         E0=0; E1=0; E2=1;
%                     end
%                 else
                    E0=exp(Vel*Re*(X0-Xmax));
                    E1=exp(Vel*Re*(X1-Xmax));
                    E2=exp(Vel*Re*(X2-Xmax));
                    DW=(E0*Y12+E1*Y20+E2*Y01);
%                 end
                k0W=(abW*(E2-E1)+acW*(E1*Y2-E2*Y1))/DW;
                k1W=(abW*(E0-E2)+acW*(E2*Y0-E0*Y2))/DW;
                k2W=(abW*(E1-E0)+acW*(E0*Y1-E1*Y0))/DW;
            end
            aW0=aW0+k0W;
            aWnb=aWnb+k1W*W(n1)+k2W*W(n2);
            aPsi0=aPsi0+(x21^2+y12^2)/Delta;
            aPsinb=aPsinb+(Psi(n1)*(y20*y12+x02*x21)+Psi(n2)*(y01*y12+x10*x21))/Delta;
            S=S+(22*W(n0)+7*W(n1)+7*W(n2))*Delta/216;
          end
        end   %Конец цикла по элементам
        if is_wall(i)==0
            Psi(i)=(-aPsinb+S)/aPsi0;
            Psi(i)=Psi(i)*QPsi+PsiOld*(1-QPsi);
            inaccPsi(i)=abs(PsiOld-Psi(i));
            W(i)=-aWnb/aW0;
            W(i)=W(i)*QW+WOld*(1-QW);
            inaccW(i)=abs(WOld-W(i));
        else
            Psi(i)=0;
            if is_wall(i)==5||is_wall(i)==6||is_wall(i)==7||is_wall(i)==8
                W(i)=0;
            else
                W(i)=-(I+aWnb)/aW0;
              if is_wall(i)==3
                  i=i;
              end
                W(i)=W(i)*QW+WOld*(1-QW);
                inaccW(i)=abs(WOld-W(i));
            end
        end
    %Конец вычисления в узле
  end   %Конец цикла по узлам 
  if sum(sum(abs(Psi)))~=0||sum(sum(abs(W)))~=0
    stp=sum(sum(inaccPsi))/sum(sum(abs(Psi)))/QPsi;
    stw=sum(sum(inaccW))/sum(sum(abs(W)))/QW;
  end
  st=max(stp,stw);
  n=n+1;
  if rem(n,50)==0
      disp(n);disp(st);
  end
end   %Конец цикла по погрешности
%Построение картинки распределения температуры
figure(2);
Tmin=min(Psi);
DeltaT=max(Psi)-Tmin; N=30; %N -количество цветов
C=zeros(1,3);
t=t';
for i=1:length(t)
    n0=t(1,i); n1=t(2,i); n2=t(3,i);
    XX=[p(1,n0) p(1,n1) p(1,n2)];
    YY=[p(2,n0) p(2,n1) p(2,n2)];
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