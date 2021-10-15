%Построение сетки в эллипсе
%Строим точки по осям
a=1.0001/1.0; 
b=1/a; 
Ni=90; Nr=120; flag=0; % a -горизонтальная ось
Rinf=2;
aMZ=0;  %0 - считать массив связей; 1 - читать с диска.
aH=1;  %0 - считать поле эллиптического проводника; 1 - читать с диска.
aTcond0=1; %0 - считать температуру и Nu теплопроводности; 1 - читать с диска.
%Находим параметр разгона по оси х
if a>b
    h=b/Ni;
else
   flag=1; h=a/Ni; a1=b; b=a; a=a1; %flag = 1 означает, что эллипс вытянут по y.
end
f=inline('q^Ni+a/h*(1-q)-1');
a2=1.00000001; b2=4; fa=f(Ni,a,h,a2); fb=f(Ni,a,h,b2);
while b2-a2>1e-7
    q=(a2+b2)/2;
    fq=f(Ni,a,h,q);
    if sign(fq)==sign(fa)
        a2=q; fa=fq;
    else
            b2=q; fb=fq;
    end
end
%Находим координаты точек на первом луче и на осях
x(1)=0; y(1)=0;
xray1(Ni)=a; yray1(Ni)=0;
for i=1:Ni
    xray1(i)=a-h*(q^(Ni-i)-1)/(q-1);
    yray1(i)=0;
    x(2+3*i*(i-1))=xray1(i); y(2+3*i*(i-1))=yray1(i);
    aell(i)=xray1(i); bell(i)=h*i;
end
%Находим координаты точек на дуге
qarc1=exp(2*log(a)/(1.5*Ni-1));
phiarc1(1)=pi/2*(qarc1-1)/(qarc1^(1.5*Ni)-1);
xarc1(1)=a; yarc1(1)=0;
for i=2:Ni*1.5
    phiarc1(i)=phiarc1(1)*(qarc1^(i)-1)/(qarc1-1);
end
for i=2:Ni*1.5+1
    xarc1(i)=a*cos(phiarc1(i-1)); yarc1(i)=b*sin(phiarc1(i-1));
    x(2+3*Ni*(Ni-1)+(i-1))=xarc1(i);
    y(2+3*Ni*(Ni-1)+(i-1))=yarc1(i);
end
%Находим координаты точек внутри первого сектора
%как точки пересечения прямых, соединяющих точки луча
%и арки с эллипсами, проведенными через точки на осях
for line=1:Ni-2
    k1(line)=(yarc1(Ni-line+1)-yray1(line))/(xarc1(Ni-line+1)-xray1(line));
    k2(line)=k1(line)*xray1(line);
    for m=1:Ni-line-1
        B=k1(line)*k2(line);
        A=k1(line)^2+(bell(m+line)/aell(m+line))^2;
        C=k2(line)^2-bell(m+line)^2;
        xline(line,m)=(B+sqrt(B^2-A*C))/A;
        yline(line,m)=bell(m+line)*sqrt(1-(xline(line,m)/aell(m+line))^2);
        x((line+m)*(line+m-1)*3+2+m)=xline(line,m);
        y((line+m)*(line+m-1)*3+2+m)=yline(line,m);
    end
end
%первый сектор заполнен
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Вычисляем координаты точек на ray2
k1ray2=yarc1(Ni+1)/xarc1(Ni+1);
for m=1:Ni
    xray2(m)=aell(m)*bell(m)/sqrt(bell(m)^2+k1ray2^2*aell(m)^2);
    yray2(m)=xray2(m)*k1ray2;
    x(2+3*m*(m-1)+m)=xray2(m); y(2+3*m*(m-1)+m)=yray2(m);
end
%Вычисляем углы точек на ark2
xarc2(1)=xarc1(Ni+1); yarc2(1)=yarc1(Ni+1);
for i=1:Ni/2
    phiarc2(i)=phiarc1(Ni+i);
end
for i=Ni/2+1:Ni-1
    phiarc2(i)=pi-phiarc2(Ni-i);
end
%вычисляем координаты точек дуги 2
for i=2:Ni
    xarc2(i)=a*cos(phiarc2(i-1)); yarc2(i)=b*sin(phiarc2(i-1));
    x(2+3*Ni*(Ni-1)+Ni+(i-1))=xarc2(i);
    y(2+3*Ni*(Ni-1)+Ni+(i-1))=yarc2(i);
end
%Вычисляем координаты внутри второго сектора
for line=1:Ni-2
    k12(line)=(yarc2(Ni-line+1)-yray2(line))/(xarc2(Ni-line+1)-xray2(line));
    k22(line)=k12(line)*xray2(line)-yray2(line);
    for m=1:Ni-line-1
        B=k12(line)*k22(line);
        A=k12(line)^2+(bell(m+line)/aell(m+line))^2;
        C=k22(line)^2-bell(m+line)^2;
        xline(line,m)=(B-sqrt(B^2-A*C))/A;
        yline(line,m)=k12(line)*xline(line,m)-k22(line);
        %yline(line,m)=bell(m+line)*sqrt(1-(xline(line,m)/aell(m+line))^2);
        n=line+m; nray2=2+3*n*(n-1)+n;
        x(nray2+m)=xline(line,m);
        y(nray2+m)=yline(line,m);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Заполняем третий сектор
for ell=1:Ni
    for m=1:ell+1
        r2=2+ell*(ell-1)*3+ell;
        r3=r2+ell;
        x(r3+m-1)=-x(r2-(m-1));
        y(r3+m-1)=y(r2-(m-1));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Заполняем нижнюю полуплоскость
for ell=1:Ni
    L=2+ell*(ell-1)*3;
    r4=L+3*ell;
    for m=1:ell*3-1
        x(r4+m)=x(r4-m);
        y(r4+m)=-y(r4-m);
    end
end
%Эллипс заполнен точками. 
%Последняя точка эллипса имеет номер 1+Ni*(Ni+1)*3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Строим массив элементов в эллипсе. Всего Ni^2 эл-тов
ell=1;
for i=1:6
    nel=i;
    t1(nel,1)=1;t1(nel,2)=i+1;
    if i~=6
        t1(nel,3)=i+2;
    else
        t1(nel,3)=2;
    end
end
nel=nel+1;
for ell=2:Ni
        for i=1:6    %находим номера точек на 6-ти лучах
            r(i)=2+(ell-2)*(ell-1)*3+(ell-1)*(i-1);
        end
        for m=1:6*(ell-1)   %цикл по точкам эллипса номер ell
            if m==1         %для первой точки два элемента
                t1(nel,1)=3*(ell-1)*(ell-2)+2+(m-1);
                t1(nel,2)=3*ell*(ell-1)+2;
                t1(nel,3)=t1(nel,2)+1;
                nel=nel+1;
                t1(nel,1)=t1(nel-1,1);
                t1(nel,2)=t1(nel-1,3);
                t1(nel,3)=t1(nel,1)+1;
                nel=nel+1;
            else            %остальные точки
                n=3*(ell-2)*(ell-1)+2+(m-1);   %номер точки
                %если точка лежит на луче (кроме первой), 3 элемента
                if n==r(2)||n==r(3)||n==r(4)||n==r(5)||n==r(6)
                    t1(nel,1)=n;
                    i=fix((m-1)/(ell-1))+1;
                    ne=2+(ell-1)*ell*3+ell*(i-1);
                    t1(nel,2)=ne-1;
                    t1(nel,3)=ne;
                    nel=nel+1;
                    t1(nel,1)=n; t1(nel,2)=ne; t1(nel,3)=ne+1;
                    nel=nel+1;
                    t1(nel,1)=n; t1(nel,2)=ne+1; 
                    if m==6*(ell-1)   %последняя точка на эллипсе ell
                        t1(nel,3)=2+(ell-2)*(ell-1)*3;
                    else
                       t1(nel,3)=n+1;
                    end
                    nel=nel+1;
                else
                    for i=1:5
                        if n>r(6)
                            dm=n-r(6); 
                            ii=6;
                        elseif (r(i)<n)&&(n<r(i+1))
                            dm=n-r(i); 
                            ii=i;
                        end
                    end
                    nrayi=2+(ell-1)*ell*3+ell*(ii-1);
                    ne=nrayi+dm;
                    t1(nel,1)=n;
                    t1(nel,2)=ne;
                    t1(nel,3)=ne+1;
                    nel=nel+1;
                    t1(nel,1)=n; t1(nel,2)=ne+1; 
                    if nel~=6*ell^2-1
                        t1(nel,3)=n+1;
                    else
                        t1(nel,3)=2+(ell-1)*(ell-2)*3;
                    end
                    nel=nel+1;
                    end
                    if m==6*(ell-1)&&nel==6*ell^2 %если точка последняя
                        t1(nel,1)=3*(ell-2)*(ell-1)+2;
                        t1(nel,2)=3*(ell-1)*ell+6*ell+1;
                        t1(nel,3)=3*(ell-1)*ell+2;
                        nel=nel+1;
                    end
            end
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%вычисляем начальные и конечные точки на радиальных линиях
%зазора
nbegin=3*Ni*(Ni-1)+2;  %первая т. на пов-ти эллипса
nend=3*Ni*(Ni+1)+1;    %посл.т. на пов-ти эллипса
%Вычисляем начальные и конечные точки лучей в зазоре, т.е на элл. и трубе
for i=1:6*Ni           %цикл вдоль пов-ти эллипса
    x0(i)=x(nbegin+i-1); y0(i)=y(nbegin+i-1);
    if i==1
        xN(i)=Rinf; yN(i)=0;
    elseif (i>1)&&(i<=1.5*Ni+1)
        xN(i)=Rinf*cos(phiarc1(i-1));
        yN(i)=Rinf*sin(phiarc1(i-1));
    elseif (i>1.5*Ni+1)&&(i<=3*Ni+1)
        xN(i)=-xN(3*Ni-i+2);
        yN(i)=yN(3*Ni-i+2);
    elseif (i>3*Ni+1)&&(i<=4.5*Ni+1)
        xN(i)=-xN(-3*Ni+i);
        yN(i)=-yN(-3*Ni+i);
    elseif (i>4.5*Ni+1)&&(i<=6*Ni+1)
        xN(i)=xN(6*Ni-i+2);
        yN(i)=-yN(6*Ni-i+2);
    end
end
%Вычисляем координаты точек в зазоре
for i=1:6*Ni
    hxi=(xN(i)-x0(i))/Nr;
    hyi=(yN(i)-y0(i))/Nr;
    for j=1:Nr
        node=nend+i+Ni*6*(j-1);
        x(node)=x0(i)+hxi*j;
        y(node)=y0(i)+hyi*j;
    end
end
%Строим массив элементов в зазоре
nel=6*Ni^2;
for i=1:Ni*6
    for j=1:Nr
        if i~=Ni*6
            nel=nel+1;
            t1(nel,1)=nbegin+(j-1)*6*Ni+(i-1);
            t1(nel,2)=nbegin+j*6*Ni+(i-1);
            t1(nel,3)=nbegin+j*6*Ni+i;
            nel=nel+1;
            t1(nel,1)=nbegin+(j-1)*6*Ni+(i-1);
            t1(nel,2)=nbegin+j*6*Ni+i;
            t1(nel,3)=nbegin+(j-1)*6*Ni+i;
        else
            nel=nel+1;
            t1(nel,1)=nend+(j-1)*6*Ni;
            t1(nel,2)=nend+j*6*Ni;
            t1(nel,3)=nbegin+j*6*Ni;
            nel=nel+1;
            t1(nel,1)=nend+(j-1)*6*Ni;
            t1(nel,2)=nbegin+j*6*Ni;
            t1(nel,3)=nbegin+(j-1)*6*Ni;
        end
    end
end
%%%%%%%%%%%%%%%%%%%
%Если a < b, меняем местами x и y.
if flag==1
    for i=1:3*Ni*(Ni+1)+1+6*Ni*Nr
        med=x(i);
        x(i)=y(i);
        y(i)=med;
    end
    a1=b; b=a; a=a1;
end
%Массив координат для основной программы
for i=1:3*Ni*(Ni+1)+1+6*Ni*Nr
    p(1,i)=x(i);
    p(2,i)=y(i);
end


% %Рисуем сетку
% figure
% axis equal
% for i=1:6*Ni^2+12*Ni*Nr
%     xxx(1)=x(t1(i,1));
%     xxx(2)=x(t1(i,2));
%     xxx(3)=x(t1(i,3));
%     xxx(4)=x(t1(i,1));
%     yyy(1)=y(t1(i,1));
%     yyy(2)=y(t1(i,2));
%     yyy(3)=y(t1(i,3));
%     yyy(4)=y(t1(i,1));
%     hold on
%     plot(xxx,yyy,'Color',[0 0 0])
% end
%Массив элементов для основной программы
t=t1';
for i=1:6*Ni^2+Nr*Ni*12
    if i<=6*Ni^2
        t(4,i)=1;  %Элемент внутри зазора
    else
        t(4,i)=2;   %Элемент в зазоре
    end
    if a<1
        ttt=t(2,i);
        t(2,i)=t(3,i);
        t(3,i)=ttt;
    end
end
% hold on
% plot(x,y,'o')
S='Mesh is constructed';
disp(S);
NombNodes=nend+6*Ni*Nr;
NombEl=6*Ni^2+Nr*Ni*12;
%p(1,i), p(2,i) - координаты x,y узлов сетки
%t(1,i),t(2,i),t(3,i) - номера узлов i-того элемента
%t(4,i) - принадлежность области (2 - зазор, 1 - внутри эллипса
%%%
%Создание массива связи: первый столюец - кол-во элтов-соседей,
%прочие - номера соседних с i-тым узлом (i - номер строки) элементов
   if aMZ==0
        MZn=zeros(NombNodes,9);
        for i=1:length(p)    %i - номер узла
            k=2;
            for j=1:length(t)    %j - номер элемента
                if i==t(1,j)||i==t(2,j)||i==t(3,j)
                    MZn(i,k)=j;
                    k=k+1;
                end
            end
            MZn(i,1)=k-2;
        end
        fmz=fopen('Aa1_0_4060.txt','wt');
        fwrite(fmz,MZn,'integer*4');
        fclose(fmz);
    else
        fmz=fopen('Aa1_0_4060.txt','rt' );
        MZn=fread(fmz,[NombNodes 9],'integer*4');
        fclose(fmz);
    end
S='Step MZn';
disp(S);
T=zeros(1,NombNodes); 
Current=25;         %Ток, А
Pr=709;
R=0.002;            %Радиус проводника, мм
Ms=44000;           %Намагниченность насыщения, А/м
Rinf=2;             %Радиус наружной трубы, безразмерный, в радиусах внутреннего
Ra=0.0; Ram=130000.0;
%%% Поле эллиптического проводника с током
% a=1.0/1.8;                  %Горизонтальная полуось
% b=1/a;                %Вертикальная полуось
if aH==0
    H=zeros(1,NombNodes);Hx=zeros(1,NombNodes); Hy=zeros(1,NombNodes);eps=zeros(1,NombNodes);
    for j=1:length(p)
        xj=p(1,j);
        yj=p(2,j);
        r=sqrt(xj^2+yj^2);
        Hx(j)=0; Hy(j)=0; eps(j)=0;
        if (xj/a)^2+(yj/b)^2>=1
            for i=1:length(t)
                x1=p(1,t(1,i));
                y1=p(2,t(1,i));
                x2=p(1,t(2,i));
                y2=p(2,t(2,i));
                x3=p(1,t(3,i));
                y3=p(2,t(3,i));
                xc=(x1+x2+x3)/3;
                yc=(y1+y2+y3)/3;
                if (xc/a)^2+(yc/b)^2<1
                    Si=((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y1-x1*y3))/2;
                    rj2=(xc-xj)^2+(yc-yj)^2;
                    Hx(j)=Hx(j)-(xj-xc)/(2*pi*rj2)*Si;
                    Hy(j)=Hy(j)+(yj-yc)/(2*pi*rj2)*Si;
                end
            end
        end
        H(j)=sqrt(Hx(j)^2+Hy(j)^2);         %Напряженность поля, без I/R
    %        eps(j)=(H(j)/(pi*a*b/(2*pi*r))-1);  Погрешность для круглого
    %        проводника
    %
    %   Обезразмеривание поля
        H(j)=H(j)*Current/(R*Ms);
    end
    fmh=fopen('H1_0_40x60.b','wt');
    fwrite(fmh,H,'double');
    fclose(fmh);
else
    fmh=fopen('H1_0_40x60.b','rt' );
    H=fread(fmh,[1 NombNodes],'double');
    fclose(fmh);
end
%Конец вычисления поля
S='Step 2 H';
disp(S);
QT=1.2; QPsi=1.2; QW=0.3;
%Начальные условия и принадлежность узла
Psi=zeros(1,NombNodes); W=zeros(1,NombNodes); T=zeros(1,NombNodes);
for i=1:length(p)
    xj=p(1,i);
    yj=p(2,i);
    r=sqrt(xj^2+yj^2);
    el=(xj/a)^2+(yj/b)^2;
    Psi(i)=0; T(i)=0;  W(i)=0;
    p(3,i)=2;   %   2 - точка принадлежит проводнику
    if r>Rinf-1e-8
        T(i)=0;
        Psi(i)=0;
        p(3,i)=0;   %   0 - точка принадлежит внешней границе зазора
    elseif (el>1-1e-8)&&(el<1+1e-8)
        T(i)=1;
        Psi(i)=0;
        p(3,i)=1;   %    1 - точка принадлежит поверхности проводника
    elseif (el>(1+1e-8))&&(r<(Rinf-1e-8))
        T(i)=0;       % начальное условие для температуры
        Psi(i)=rand(1)*0.001;               %начальное условие для функции тока
        W(i)=rand(1)*0.001;                     %начальное условие для вихря
        p(3,i)=3;   %    3 - точка принадлежит зазору
    end
end
%Дополняем массив t восприимчивостью
for i=1:length(t)      %элементы, для которых t(4,i)==2 находятся в зазоре
    if t(4,i)==2       %пятую строку массива t заполняем значениями
        t(5,i)=5;      %начальной восприимчивости жидкости: chi=5 - в зазоре
    else
        t(5,i)=0;      %   chi=0 - в проводнике
    end
end
S='Step 3, Initial conditions';
disp(S);
%Вычисление температуры в режиме теплопроводности
if aTcond0==0

l=0;
disp(l);
inaccuracy=10; inacc=0; inaccT=zeros(1,NombNodes); 
inaccPsi=zeros(1,NombNodes); inaccW=zeros(1,NombNodes); inaccF=zeros(1,NombNodes);
Tsum=0; 
inaccuracy=10; lw=0;
while(inaccuracy>1e-7)  %Цикл T теплопроводность
   for i=1:length(p)       % i цикл по узлам 
       if p(3,i)==3        % 3 - внутри зазора
            a0T=0; anbT=0;
            for j=1:MZn(i,1)    % цикл по эл-там вокруг узла
                nel=MZn(i,j+1);  %номер соседнего эл-та
                n0=t(1,nel);n1=t(2,nel);n2=t(3,nel);
                if i==n1
                    n1=n2;n2=n0;n0=i;
                end
                if i==n2
                    n2=n1;n1=n0;n0=i;
                end
                x0=p(1,n0);y0=p(2,n0);
                x1=p(1,n1);y1=p(2,n1);
                x2=p(1,n2);y2=p(2,n2);
                x01=x0-x1;y01=y0-y1;
                x12=x1-x2;y12=y1-y2;
                x20=x2-x0;y20=y2-y0;
                Dpsi=(x1*y2-x2*y1)+(x2*y0-x0*y2)+(x0*y1-x1*y0);
                a0T=a0T-0.5/Dpsi*(y12*y12+x12*x12);
                anbT=anbT-0.5/Dpsi*(T(n1)*(y12*y20+x12*x20)+T(n2)*(y12*y01+x12*x01));
            end   %Конец цикла по элементам вокруг узла
            TOld=T(i);
            T(i)=-anbT/a0T*QT+TOld*(1-QT);
            inaccT(i)=abs((TOld-T(i))/T(i));
       end    %Конец if по внутренним точкам
   end        %Конец цикла по узлам
   if sum(abs(T))==0
       inaccTS=0;
   else
       inaccTS=sum(inaccT)/sum(abs(T));
   end
   lw=lw+1;
   inaccuracy=min(inaccuracy,inaccTS);
   if rem(lw,100)==0        %Вывод на экран промежуточных данных о сходимости
       disp(lw);
       disp(inaccTS);
   end
end     %Конец цикла T теплопроводность
disp(inaccuracy);
        fmTu=fopen('Tcond01_0_40x60.txt','wt');
        fwrite(fmTu,T,'double');
        fclose(fmTu);
else
        fmTu=fopen('Tcond01_0_40x60.txt','rt' );
        T=fread(fmTu,[NombNodes 1],'double');
        fclose(fmTu);
end


%Вычисление интегрального Нуссельта       
Nu=0;
for i=1:length(p)
    if p(3,i)==1
        for j=1:MZn(i,1)
          nel=MZn(i,j+1);  %номер соседнего эл-та
          n0=t(1,nel);n1=t(2,nel);n2=t(3,nel);
          if i==n1
              n1=n2;n2=n0;n0=i;
          end
          if i==n2
              n2=n1;n1=n0;n0=i;
          end
          if (p(3,n2)==1)&&(t(4,nel)==2)
                x0=p(1,n0);y0=p(2,n0);
                x1=p(1,n1);y1=p(2,n1);
                x2=p(1,n2);y2=p(2,n2);
                x01=x0-x1;y01=y0-y1;
                x12=x1-x2;y12=y1-y2;
                x20=x2-x0;y20=y2-y0;
                Dpsi=(x1*y2-x2*y1)+(x2*y0-x0*y2)+(x0*y1-x1*y0);
                BT=-(T(n0)*x12+T(n1)*x20+T(n2)*x01)/Dpsi;
                AT=(T(n0)*y12+T(n1)*y20+T(n2)*y01)/Dpsi;
                s=sqrt(x20^2+y20^2);
                nx=(-y20)/s;
                ny=(x20)/s;
                Nuloc=(AT*nx+BT*ny);
                Nu=Nu+Nuloc*s;
          end
        end
    end
end
SS='Nu0 =';
disp(SS)
disp(Nu);          %Вывод на экран Нуссельта для теплопроводности

l=0; disp(l);
inaccuracy=10; inacc=zeros(1,NombNodes); inaccT=zeros(1,NombNodes); 
inaccPsi=zeros(1,NombNodes); inaccW=zeros(1,NombNodes); %inaccF=0; H=0;
%Вычисление T, Psi, W
% Внешний цикл по T, Psi, W
Tsum=0; Psum=0; Wsum=0; ltpw=1; acc=10; 
%        
inaccuracy=10; lw=0;
while(inaccuracy>1e-8)  %Цикл T,P,W. Начало цикла по ф-ии тока - вихрю и т-ре
   for i=1:length(p)       % i цикл по узлам 
     inaccW(i)=0; inaccT(i)=0;   inaccPsi(i)=0;
     if p(3,i)==0||p(3,i)==1||p(3,i)==3     %условие принадлежности точки "зазору"
            a0W=0; anbW=0; a0T=0; anbT=0; 
            anbPsi=0;  a0Psi=0; 
            for j=1:MZn(i,1)    % цикл по эл-там вокруг узла
                nel=MZn(i,j+1);  %номер соседнего эл-та
                if t(4,nel)==2        % элемент внутри зазора, не в проводнике
                n0=t(1,nel);n1=t(2,nel);n2=t(3,nel);
                if i==n1
                    n1=n2;n2=n0;n0=i;
                end
                if i==n2
                    n2=n1;n1=n0;n0=i;
                end
                x0=p(1,n0);y0=p(2,n0);
                x1=p(1,n1);y1=p(2,n1);
                x2=p(1,n2);y2=p(2,n2);
                x01=x0-x1;y01=y0-y1;
                x12=x1-x2;y12=y1-y2;
                x20=x2-x0;y20=y2-y0;
                Dpsi=(x1*y2-x2*y1)+(x2*y0-x0*y2)+(x0*y1-x1*y0);
                DAPsi=Psi(n0)*y12+Psi(n1)*y20+Psi(n2)*y01;
                DBPsi=-(Psi(n0)*x12+Psi(n1)*x20+Psi(n2)*x01);
                APsi=DAPsi/Dpsi; BPsi=DBPsi/Dpsi;
                Uel=sqrt(APsi*APsi+BPsi*BPsi);
                if abs(Uel)<3e-15
                    Uel=1e-8; sinalpha=0; cosalpha=1;
                end
                xc=(x0+x1+x2)/3; yc=(y0+y1+y2)/3;
                sinalpha=-APsi/Uel; cosalpha=BPsi/Uel;    
                X0=(x0-xc)*cosalpha+(y0-yc)*sinalpha; 
                X1=(x1-xc)*cosalpha+(y1-yc)*sinalpha;  
                X2=(x2-xc)*cosalpha+(y2-yc)*sinalpha; 
                Y0=-(x0-xc)*sinalpha+(y0-yc)*cosalpha;
                Y1=-(x1-xc)*sinalpha+(y1-yc)*cosalpha; 
                Y2=-(x2-xc)*sinalpha+(y2-yc)*cosalpha;
                X01=X0-X1;Y01=Y0-Y1;
                X12=X1-X2;Y12=Y1-Y2;
                X20=X2-X0;Y20=Y2-Y0;
                Xmax=max(X0,X1); Xmax=max(Xmax,X2);
                Ymax=max(Y0,Y1); Ymax=max(Ymax,Y2);
                if p(3,i)==3        %узел внутри зазора
                    anbPsi=anbPsi+1/Dpsi*(Psi(n1)*(y12*y20+x12*x20)+Psi(n2)*(y12*y01+x12*x01))-Dpsi/108*(W(n0)*22+W(n1)*7+W(n2)*7);
                    a0Psi=a0Psi+1/Dpsi*(y12*y12+x12*x12);
                    E0W=exp(Uel*(X0-Xmax)/Pr);
                    E1W=exp(Uel*(X1-Xmax)/Pr);
                    E2W=exp(Uel*(X2-Xmax)/Pr);
                    DW=E0W*Y12+E1W*Y20+E2W*Y01;
                    aBW=Uel*0.125*Y12*(Y1+2*Y0+Y2)/Pr+0.5*X12;
                    aCW=Uel*0.5*Y12/Pr;
                    k1W=aBW*(E0W-E2W)+aCW*(E2W*Y0-E0W*Y2);
                    k2W=aBW*(E1W-E0W)+aCW*(E0W*Y1-E1W*Y0);
                    k0t=(y2-y1)/24;
                    k1t=-(15*y0+9*y1+12*y2)/72;
                    k2t=(15*y0+12*y1+9*y2)/72;
                    if lw>1
                        k0Ht=(H(n2)-H(n1))/24;
                        k1Ht=-(15*H(n0)+9*H(n1)+12*H(n2))/72;
                        k2Ht=(15*H(n0)+12*H(n1)+9*H(n2))/72;
                    end
                    anbW=anbW+((W(n1)*k1W+W(n2)*k2W)/DW-Ra*(k0t*T(n0)+k1t*T(n1)+k2t*T(n2)));
                    if lw>1
                        anbW=anbW+Ram*t(5,nel)*H(i)/(1+t(5,nel)*H(i))*(k0Ht*T(n0)+k1Ht*T(n1)+k2Ht*T(n2));
                    end
                    a0W=a0W+(aBW*(E2W-E1W)+aCW*(E1W*Y2-E2W*Y1))/DW;
                    E0T=exp(Uel*(X0-Xmax));
                    E1T=exp(Uel*(X1-Xmax));
                    E2T=exp(Uel*(X2-Xmax));
                    DT=E0T*Y12+E1T*Y20+E2T*Y01;
                    aBT=Uel*0.125*Y12*Y0+0.5*X12;
                    aCT=Uel*0.5*Y12;
                    a0T=a0T+(aBT*(E2T-E1T)+aCT*(E1T*Y2-E2T*Y1))/DT;
                    k1T=aBT*(E0T-E2T)+aCT*(E2T*Y0-E0T*Y2);
                    k2T=aBT*(E1T-E0T)+aCT*(E0T*Y1-E1T*Y0);
                    anbT=anbT+(T(n1)*k1T+T(n2)*k2T)/DT;
                elseif p(3,i)==0||p(3,i)==1       %если узел лежит на одной из границ зазора
                    if p(3,i)==0                %наружная граница
                        x1loc=x01*y0/Rinf-y01*x0/Rinf;
                        x2loc=-x20*y0/Rinf+y20*x0/Rinf;
                        y1loc=x01*x0/Rinf+y01*y0/Rinf;
                        y2loc=-x20*x0/Rinf-y20*y0/Rinf;
                    end
                    if p(3,i)==1                %внутренняя граница
                        nx=b^2*x0/sqrt(b^4*x0^2+a^4*y0^2);
                        ny=a^2*y0/sqrt(b^4*x0^2+a^4*y0^2);
                        x1loc=-x01*ny+y01*nx;
                        y1loc=-x01*nx-y01*ny;
                        x2loc=x20*ny-y20*nx;
                        y2loc=x20*nx+y20*ny;
                    end
                    xc=(x1loc+x2loc)/3; yc=(y1loc+y2loc)/3;
                    xa=0.5*x1loc; ya=0.5*y1loc;
                    xb=0.5*x2loc; yb=0.5*y2loc;
                    Ia=(yb*yb-ya*ya)-(xb*xb-xa*xa);
                    Ib=(ya+yc)*(xc-xa)+(yc+yb)*(xb-xc);
                    D=0.5*(y1loc*y1loc*x2loc*y2loc-y2loc*y2loc*x1loc*y1loc);
                    Db=Psi(n1)*x2loc*y2loc-Psi(n2)*x1loc*y1loc;
                    Da=0.5*(Psi(n2)*y1loc*y1loc-Psi(n1)*y2loc*y2loc);
                    Di=x1loc*y2loc-x2loc*y1loc;
                    anbW=anbW+(Da*Ia-Db*Ib)/D+(7*W(n1)+7*W(n2))*Di/108;
                    a0W=a0W+11/54*Di;
                 end             % if узел внутри цилиндра
              end            %конец проверки принадлежности эл-та зазору
            end            % Конец цикла по элементам вокруг узла
            PsiOld=Psi(i); WOld=W(i); TOld=T(i);
            if p(3,i)==0||p(3,i)==1
                if p(3,i)==0
                    T(i)=0;
                elseif p(3,i)==1
                    T(i)=1;
                end
                Psi(i)=0;
                 W(i)=-anbW/a0W*QW+WOld*(1-QW);             
            else
                if p(3,i)==3
                    T(i)=-anbT/a0T*QT+TOld*(1-QT);
                    inaccT(i)=abs(TOld-T(i));
                    Psi(i)=-anbPsi/a0Psi*QPsi+PsiOld*(1-QPsi);
                    inaccPsi(i)=abs(PsiOld-Psi(i));
                    W(i)=-anbW/a0W*QW+WOld*(1-QW);
                end
            end
            inaccW(i)=abs(WOld-W(i));
     end             %Конец условия "зазор"
   end          %конец цикла по всем узлам
        if sum(abs(Psi))==0
            inaccPsiS=0;
        else
            inaccPsiS=sum(inaccPsi)/sum(abs(Psi));
        end
        if sum(abs(W))==0
            inaccWS=0;
        else
            inaccWS=sum(inaccW)/sum(abs(W));
        end
        if sum(abs(T))==0
            inaccTS=0;
        else
            inaccTS=sum(inaccT)/sum(abs(T));
        end
    in2=max(inaccPsiS,inaccWS);
    inaccuracy=max(inaccTS,in2);
    lw=lw+1;
    if rem(lw,50)==0
        disp(inaccPsiS);
        disp(inaccWS);
        disp(inaccTS);
        disp(inaccuracy);
        disp(lw);
    end
end     %Конец цикла T,P,W по ф-ии тока - вихрю и т-ре
%Вычисление интегрального Нуссельта       
Nu1=0;
for i=1:length(p)
    if p(3,i)==1
        for j=1:MZn(i,1)
          nel=MZn(i,j+1);  %номер соседнего эл-та
          n0=t(1,nel);n1=t(2,nel);n2=t(3,nel);
          if i==n1
              n1=n2;n2=n0;n0=i;
          end
          if i==n2
              n2=n1;n1=n0;n0=i;
          end
          if (p(3,n2)==1)&&(t(4,nel)==2)
                x0=p(1,n0);y0=p(2,n0);
                x1=p(1,n1);y1=p(2,n1);
                x2=p(1,n2);y2=p(2,n2);
                x01=x0-x1;y01=y0-y1;
                x12=x1-x2;y12=y1-y2;
                x20=x2-x0;y20=y2-y0;
                Dpsi=(x1*y2-x2*y1)+(x2*y0-x0*y2)+(x0*y1-x1*y0);
                BT=-(T(n0)*x12+T(n1)*x20+T(n2)*x01)/Dpsi;
                AT=(T(n0)*y12+T(n1)*y20+T(n2)*y01)/Dpsi;
                s=sqrt(x20^2+y20^2);
                nx=(-y20)/s;
                ny=(x20)/s;
                Nuloc=(AT*nx+BT*ny);
                Nu1=Nu1+Nuloc*s;
          end
        end
    end
end
disp(Nu1);          %Вывод на экран Нуссельта
%Картинка линий тока. Если вместо Psi вставить Т, будут изотермы
% t1=[t(1,:);t(2,:);t(3,:);t(4,:)];
% p1=[p(1,:);p(2,:)];
% pdeplot(p1,e,t1,'xydata',Psi,'colorbar','on','contour','on','levels',20)
    