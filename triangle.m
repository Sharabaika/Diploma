%Построение сетки с треугольными элементами в прямоугольнике
%Построение массива координат p, массива элементов t модулем Circles.m
global p t;
aMZ=0;     %Флаг, 0 - считаем, 1 - читаем массив MZn
H0 = 1;     %Внешнее поле в единицах Ms
chi0=3;     %Начальная восприимчивость
QF=1.6;     %Параметр релаксации
%%%%%%%%%%%%%%%%%%%
NumbNodes=length(p);
NumbElements=length(t);
%Рисуем сетку
figure(1)
axis equal
hold on
t=t';       %транспонируем, т.к. модуль создает массив t горизонтальным
for i=1:NumbElements
    xxx(1)=p(1,t(i,1));
    xxx(2)=p(1,t(i,2));
    xxx(3)=p(1,t(i,3));
    xxx(4)=p(1,t(i,1));
    yyy(1)=p(2,t(i,1));
    yyy(2)=p(2,t(i,2));
    yyy(3)=p(2,t(i,3));
    yyy(4)=p(2,t(i,1));
    plot(xxx,yyy,'Color',[0 0 0])
end

%Создание массива связи: первый столюец - кол-во элтов-соседей,
%прочие - номера соседних с i-тым узлом (i - номер строки) элементов
NumbNodes=length(p);   
if aMZ==0
        MZn=zeros(NumbNodes,9);
        for i=1:NumbNodes    %i - номер узла
            k=2;
            for j=1:NumbElements    %j - номер элемента
                if i==t(j,1)||i==t(j,2)||i==t(j,3)
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
S='Step MZn'



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Вычисление потенцила поля
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Заполнение массива граничных точек
is_wall=zeros(1,NumbNodes);
F=zeros(1,NumbNodes); %Массив потенциала
for i=1:NumbNodes
    Ri=sqrt(p(1,i)^2+p(2,i)^2);
    F(i)=H0*p(2,i);     %Задаем начальное распределение потенциала
    if abs(Ri-10)<1e-8
        is_wall(i)=1;  %Узлы на внешней границе
    else
        is_wall(i)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Начальное распределение магнитной проницаемости по элементам
mu=zeros(1,NumbElements);
H=zeros(1,NumbElements); %Массив напряженности
for i=1:NumbElements
    H(i)=H0;            %Начальное распределение поле однородное
    if t(i,4)==2
        mu(i)=1;
    elseif t(i,4)==1
        mu(i)=1+chi0;  %/(1+chi0*H(i));
    end
end
inaccF=zeros(1,NumbNodes); 
inaccuracy=10; lw=0;
while(inaccuracy>1e-7)  %Цикл итераций
   for i=1:NumbNodes       % i цикл по узлам 
            a0F=0; anbF=0;
            FOld=F(i);
            if is_wall(i)==0
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
                Delta=x10*y20-x02*y01;
                a0F=a0F+mu(nel)*0.5/Delta*(y12*y12+x21*x21);
                anbF=anbF+mu(nel)*0.5/Delta*(F(n1)*(y12*y20+x21*x02)+F(n2)*(y12*y01+x21*x10));
              end   %Конец цикла по элементам вокруг узла
              F(i)=-anbF/a0F*QF+FOld*(1-QF);
              inaccF(i)=abs(FOld-F(i));
            elseif is_wall(i)==1
                F(i)=H0*p(2,i);
            end    
   end        %Конец цикла по узлам
%Считаем поле и проницаемость на каждом элементе
    for i=1:NumbElements
        n0=t(i,1);n1=t(i,2);n2=t(i,3);
        x0=p(1,n0);y0=p(2,n0);
        x1=p(1,n1);y1=p(2,n1);
        x2=p(1,n2);y2=p(2,n2);
        x10=x1-x0;y01=y0-y1;
        x21=x2-x1;y12=y1-y2;
        x02=x0-x2;y20=y2-y0;
        Delta=x10*y20-x02*y01;
        DeltaA=F(n0)*y12+F(n1)*y20+F(n2)*y01;
        DeltaB=F(n0)*x21+F(n1)*x02+F(n2)*x10;
        Hx=DeltaA/Delta;
        Hy=DeltaB/Delta;
        H(i)=sqrt(Hx^2+Hy^2);
        if t(i,4)==2
            mu(i)=1;   %проницаемость вне круга радиусом 1
        elseif t(i,4)==1
            mu(i)=1+chi0; %/(1+chi0*H(i)); %проницаемость внутри круга радиусом 1
        end
    end
%Конец вычисления поля и проницаемости на ЭЛЕМЕНТАХ (не в узлах!) 
    if sum(abs(F))==0
       inaccuracy=1;
   else
       inaccuracy=sum(inaccF)/sum(abs(F));
   end
   lw=lw+1;
   if rem(lw,100)==0        %Вывод на экран промежуточных данных о сходимости
       disp(lw);
       disp(inaccuracy);
   end
end     %Конец итерационного цикла
figure(2);
DeltaT=max(H)-min(H); 
N=30; %N -количество цветов
C=zeros(1,3);
t=t';
for i=1:NumbElements
    n0=t(1,i); n1=t(2,i); n2=t(3,i);
    XX=[p(1,n0) p(1,n1) p(1,n2)];
    YY=[p(2,n0) p(2,n1) p(2,n2)];
    TT=(H(i)-min(H))/DeltaT;
    %h=1/(N/3);
    if TT<1/3
        C(1)=0; C(2)=3*TT; C(3)=1;
    elseif (TT>=1/3)&&(TT<2/3)
        C(1)=3*(TT-1/3); C(2)=1; C(3)=1-3*(TT-1/3);
        %C(1)=0; C(2)=1-3*(TT-1/3); C(3)=0;
    elseif (TT>=2/3)&&(TT<=1)
        C(1)=1; C(2)=1-3*(TT-2/3); C(3)=0;
        %C(1)=-9*(TT-1)^2+1; C(2)=9*(TT-1)^2; C(3)=0;
        if (TT<2*DeltaT/3+0.1*DeltaT/3)
            i=i;
        end
    end
    fill(XX,YY,C,'LineStyle','none');
    hold on;
end