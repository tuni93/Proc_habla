%Guia 6
%Funes Pablo Nicolas
%Padron 94894

clear all
close all

%Lectura de los archivos

Muestras_A=load('a.txt');
Muestras_O=load('o.txt');
Muestras_U=load('u.txt');

Muestras_A_desordenadas(:,1)=Muestras_A(randperm(length(Muestras_A(:,1))));
Muestras_A_desordenadas(:,2)=Muestras_A(randperm(length(Muestras_A(:,2))));

Muestras_O_desordenadas(:,1)=Muestras_O(randperm(length(Muestras_O(:,1))));
Muestras_O_desordenadas(:,2)=Muestras_O(randperm(length(Muestras_O(:,2))));

Muestras_U_desordenadas(:,1)=Muestras_U(randperm(length(Muestras_U(:,1))));
Muestras_U_desordenadas(:,2)=Muestras_U(randperm(length(Muestras_U(:,2))));

N_max=25;%Nro maximo de posibles muestras de una clase a tomar
Nro_muestras_aleatorias_A=round(rand()*(N_max-1)+1)
Nro_muestras_aleatorias_O=round(rand()*(N_max-1)+1)
Nro_muestras_aleatorias_U=round(rand()*(N_max-1)+1)
Nro_total_muestras=Nro_muestras_aleatorias_A+Nro_muestras_aleatorias_O+Nro_muestras_aleatorias_U;
Muestras_aleatorias=[Muestras_A_desordenadas([1:Nro_muestras_aleatorias_A],:);Muestras_O_desordenadas([1:Nro_muestras_aleatorias_O],:);Muestras_U_desordenadas([1:Nro_muestras_aleatorias_U],:)];
Vector_unos=ones(1,Nro_total_muestras);
Punto_medio=(Vector_unos*Muestras_aleatorias)/Nro_total_muestras;

figure(1)
hold on
plot(Muestras_aleatorias(:,1),Muestras_aleatorias(:,2),'ro');
plot(Punto_medio(1,1)*Vector_unos,Punto_medio(1,2)*Vector_unos,'bo');
%%%%%%%%%%%%%%%%%
%RECTA ALEATORIA%
%%%%%%%%%%%%%%%%%

pendiente_1=rand()*2*pi;
angulo=120*pi/180;
X_recta_1=[0:500/Nro_total_muestras:500-1];
Y_recta_1=X_recta_1*tan(pendiente_1);
matriz_rotacion=[cos(angulo) -sin(angulo);sin(angulo) cos(angulo)];
X_recta_2 =matriz_rotacion(1,:)*[X_recta_1;Y_recta_1];
Y_recta_2 =matriz_rotacion(2,:)*[X_recta_1;Y_recta_1];
X_recta_3 =matriz_rotacion(1,:)*[X_recta_2;Y_recta_2];
Y_recta_3 =matriz_rotacion(2,:)*[X_recta_2;Y_recta_2];


hold on
plot(X_recta_1+Punto_medio(1,1),Y_recta_1+Punto_medio(1,2));
plot(X_recta_2+Punto_medio(1,1),Y_recta_2+Punto_medio(1,2));
plot(X_recta_3+Punto_medio(1,1),Y_recta_3+Punto_medio(1,2));


%Puntos Iniciales en cada area aleatorios

X_punto_1=100*cos(pendiente_1/2);
Y_punto_1=100*sin(pendiente_1/2);
X_punto_2 =matriz_rotacion(1,:)*[X_punto_1;Y_punto_1];
Y_punto_2 =matriz_rotacion(2,:)*[X_punto_1;Y_punto_1];
X_punto_3 =matriz_rotacion(1,:)*[X_punto_2;Y_punto_2];
Y_punto_3 =matriz_rotacion(2,:)*[X_punto_2;Y_punto_2];

X_punto_1=X_punto_1+Punto_medio(1,1);
Y_punto_1=Y_punto_1+Punto_medio(1,2);
X_punto_2=X_punto_2+Punto_medio(1,1);
Y_punto_2=Y_punto_2+Punto_medio(1,2);
X_punto_3=X_punto_3+Punto_medio(1,1);
Y_punto_3=Y_punto_3+Punto_medio(1,2);

hold on
plot(X_punto_1*Vector_unos,Y_punto_1*Vector_unos,'gx');
plot(X_punto_2*Vector_unos,Y_punto_2*Vector_unos,'c*');
plot(X_punto_3*Vector_unos,Y_punto_3*Vector_unos,'m+');
legend('Muestras aleatorias','MEDIA GRAL','RECTA 1','RECTA 2','RECTA 3','Media clase 1','Media clase 2','Media clase 3');


Media_A=[X_punto_1 Y_punto_1];
Media_O=[X_punto_2 Y_punto_2];
Media_U=[X_punto_3 Y_punto_3];


Nro_muestras_A=40;
Nro_muestras_O=40;
Nro_muestras_U=40;



%Armo un vector con todas las muestras
Muestras_totales=[Muestras_A_desordenadas([1:Nro_muestras_A],:);Muestras_O_desordenadas([1:Nro_muestras_O],:);Muestras_U_desordenadas([1:Nro_muestras_U],:)];

%Calculo la cantidad de muestras
Nro_total_muestras=length(Muestras_totales);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%K-MEANS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Armo un vector que tenga la clase de cada muestra
clase_cercana=zeros(1,Nro_total_muestras);
%Condiciones para la primera iteracion
distorsion_paso_previo=0;
distorsion_antigua=0;
distorsion_actualizada=1;
j=1;
k=2;
%criterio de iteracion:itero hasta que la distorsion varie menos de un 10%
while(abs((distorsion_actualizada-distorsion_antigua)/distorsion_actualizada)>0.001) && (distorsion_paso_previo~=distorsion_actualizada)
    distorsion_paso_previo=distorsion_antigua;
    distorsion_antigua=distorsion_actualizada;
    distorsion_actualizada=0;
    i=1;
    cantidad_clase_1=0;
    cantidad_clase_2=0;
    cantidad_clase_3=0;
    Media_actualizada_A=[0 0];
    Media_actualizada_O=[0 0];
    Media_actualizada_U=[0 0];
    Catalogadas_A=zeros(100,2);
    Catalogadas_O=zeros(100,2);
    Catalogadas_U=zeros(100,2);
    while(i<=Nro_total_muestras)
        %clasifico
        [aux,distancia]=Clasificacion_clase(Muestras_totales(i,:),Media_A,Media_O,Media_U);
        clase_cercana(1,i)=aux;
        distorsion_actualizada=distorsion_actualizada+distancia;
        %Cuento la cantidad de cada clase y hago la suma de todos los elementos
        %de una clase
        if(clase_cercana(1,i)==1)
            Media_actualizada_A=Media_actualizada_A+Muestras_totales(i,:);
            cantidad_clase_1=cantidad_clase_1+1;
            Catalogadas_A(cantidad_clase_1,:)=Muestras_totales(i,:)
        elseif(clase_cercana(1,i)==2)
            Media_actualizada_O=Media_actualizada_O+Muestras_totales(i,:);
            cantidad_clase_2=cantidad_clase_2+1;
            Catalogadas_O(cantidad_clase_2,:)=Muestras_totales(i,:)
        elseif(clase_cercana(1,i)==3)
            Media_actualizada_U=Media_actualizada_U+Muestras_totales(i,:);
            cantidad_clase_3=cantidad_clase_3+1;
            Catalogadas_U(cantidad_clase_3,:)=Muestras_totales(i,:)
        end
        i=i+1;  
    end

%Dividiendo la suma por la cantidad de cada clase obtengo las nuevas
%medias y divido la distorsion en la cantidad de muestras obtengo la
%distorsion.

Media_actualizada_A=Media_actualizada_A/cantidad_clase_1;
Media_actualizada_O=Media_actualizada_O/cantidad_clase_2;
Media_actualizada_U=Media_actualizada_U/cantidad_clase_3;

Media_A=Media_actualizada_A;
Media_O=Media_actualizada_O;
Media_U=Media_actualizada_U;


distorsion_actualizada=distorsion_actualizada/Nro_total_muestras;
distorsion_vector(1,j)=distorsion_antigua;
distorsion_vector(2,j)=distorsion_actualizada;
j=j+1;
%Agrego al grafico la media
figure(k)
hold on
plot(Vector_unos*Media_A(1,1),Vector_unos*Media_A(1,2),'kx');
plot(Vector_unos*Media_O(1,1),Vector_unos*Media_O(1,2),'k*');
plot(Vector_unos*Media_U(1,1),Vector_unos*Media_U(1,2),'k+');
plot(Catalogadas_A([1:cantidad_clase_1],1),Catalogadas_A([1:cantidad_clase_1],2),'ro')
plot(Catalogadas_O([1:cantidad_clase_2],1),Catalogadas_O([1:cantidad_clase_2],2),'bo')
plot(Catalogadas_U([1:cantidad_clase_3],1),Catalogadas_U([1:cantidad_clase_3],2),'go')
k=k+1;

end


%Grafico la distorsion%
%figure(2)
%hold on
%plot(distorsion_vector(2,:),'or')
%legend('Distorsion actual');


%Calculo los pesos de cada clase
pi_A=cantidad_clase_1/Nro_total_muestras;
pi_O=cantidad_clase_2/Nro_total_muestras;
pi_U=cantidad_clase_3/Nro_total_muestras;

%Calculo las varianzas
%Separo las muestras y las guardo en vectores.
i=1;
k1=1;
k2=1;
k3=1;
while(i<=Nro_total_muestras)
   if(clase_cercana(1,i)==1)
       Muestras_cl_A(k1,:)=Muestras_totales(i,:);
       k1=k1+1;
   elseif(clase_cercana(1,i)==2)
       Muestras_cl_O(k2,:)=Muestras_totales(i,:);
       k2=k2+1;
   elseif(clase_cercana(1,i)==3)
       Muestras_cl_U(k3,:)=Muestras_totales(i,:);
       k3=k3+1;
   end
   i=i+1;
end

Matriz_unos_columnas_A=ones(length(Muestras_cl_A),1);
Matriz_unos_columnas_O=ones(length(Muestras_cl_O),1);
Matriz_unos_columnas_U=ones(length(Muestras_cl_U),1);

Vector_media_A=Matriz_unos_columnas_A*Media_A;
Vector_media_O=Matriz_unos_columnas_O*Media_O;
Vector_media_U=Matriz_unos_columnas_U*Media_U;

Varianza_A=(transpose(Muestras_cl_A-Vector_media_A))*(Muestras_cl_A-Vector_media_A)/length(Muestras_cl_A);
Varianza_O=(transpose(Muestras_cl_O-Vector_media_O))*(Muestras_cl_O-Vector_media_O)/length(Muestras_cl_O);
Varianza_U=(transpose(Muestras_cl_U-Vector_media_U))*(Muestras_cl_U-Vector_media_U)/length(Muestras_cl_U);

%Varianza=(cantidad_clase_1*Varianza_A+cantidad_clase_2*Varianza_O+cantidad_clase_3*Varianza_U)/(cantidad_clase_1+cantidad_clase_2+cantidad_clase_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%PRUEBA%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Muestras_evaluacion=[Muestras_A_desordenadas([Nro_muestras_A+1:end],:);Muestras_O_desordenadas([Nro_muestras_O+1:end],:);Muestras_U_desordenadas([Nro_muestras_U+1:end],:)];
Longitud_test=length(Muestras_evaluacion);
i=1;
while(i<=Longitud_test)
    i;%Agrego un 1 en log pk
    GKA=(-0.5)*log(det(Varianza_A))-0.5*(Muestras_evaluacion(i,:)-Media_A)*(Varianza_A^-1)*transpose(Muestras_evaluacion(i,:)-Media_A)+log(pi_A);
    GKO=(-0.5)*log(det(Varianza_O))-0.5*(Muestras_evaluacion(i,:)-Media_O)*(Varianza_O^-1)*transpose(Muestras_evaluacion(i,:)-Media_O)+log(pi_O);
    GKU=(-0.5)*log(det(Varianza_U))-0.5*(Muestras_evaluacion(i,:)-Media_U)*(Varianza_U^-1)*transpose(Muestras_evaluacion(i,:)-Media_U)+log(pi_U);

   valores=[GKA GKO GKU];
   respuesta(i) = find(max(valores) == valores);
   i=i+1;
end
Matriz=ones(1,length(Muestras_evaluacion));

%figure(3)
%hold on
%plot(Muestras_evaluacion(:,1),Muestras_evaluacion(:,2),'ro')
%plot(Matriz*Media_A(1,1),Matriz*Media_A(1,2),'rx');
%plot(Matriz*Media_O(1,1),Matriz*Media_O(1,2),'g*');
%plot(Matriz*Media_U(1,1),Matriz*Media_U(1,2),'b+');
