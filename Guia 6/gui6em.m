%Algoritmo EM
%Funes Pablo Nicolas
%Padron 94894
close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%USO EL MISMO CODIGO QUE KMEANS CON BOOTSTRAP%%%%%%%%%%%%%%%
%%%%%%%%%%PARA DARLE CONDICIONES INICIALES%%%%%%%%%%%%%%%%%%%%%%%%%%%

Muestras_A=load('a.txt');
Muestras_O=load('o.txt');
Muestras_U=load('u.txt');

%Defino el numero de muestras evaluacion
Nro_muestras_evaluacion=10;
Nro_muestras_entrenamiento=5;
%Calculos la cantidad de muestras de entrenamiento

Nro_muestras_A=length(Muestras_A)-Nro_muestras_evaluacion;
Nro_muestras_O=length(Muestras_O)-Nro_muestras_evaluacion;
Nro_muestras_U=length(Muestras_U)-Nro_muestras_evaluacion;

%Desordeno las muestras para simular la aleatoriedad
F1_A=Muestras_A(:,1);
F2_A=Muestras_A(:,2);
F1_O=Muestras_O(:,1);
F2_O=Muestras_O(:,2);
F1_U=Muestras_U(:,1);
F2_U=Muestras_U(:,2);

Muestras_A_desordenadas(:,1)=F1_A(randperm(length(F1_A)));
Muestras_A_desordenadas(:,2)=F2_A(randperm(length(F2_A)));

Muestras_O_desordenadas(:,1)=F1_O(randperm(length(F1_O)));
Muestras_O_desordenadas(:,2)=F2_O(randperm(length(F2_O)));

Muestras_U_desordenadas(:,1)=F1_U(randperm(length(F1_U)));
Muestras_U_desordenadas(:,2)=F2_U(randperm(length(F2_U)));

%Tomo parte de las muestras para entrenamiento conociendo su respectiva
%clase
%Considero como muestras de entrenamiento las primeras
%Nro_muestras_entrenamiento
Matriz_unos_entrenamiento=ones(1,Nro_muestras_entrenamiento);

%Calculo la media de las muestras de entrenamiento
Media_A=(Matriz_unos_entrenamiento*Muestras_A_desordenadas([1:Nro_muestras_entrenamiento],:))/Nro_muestras_entrenamiento;
Media_O=(Matriz_unos_entrenamiento*Muestras_O_desordenadas([1:Nro_muestras_entrenamiento],:))/Nro_muestras_entrenamiento;
Media_U=(Matriz_unos_entrenamiento*Muestras_U_desordenadas([1:Nro_muestras_entrenamiento],:))/Nro_muestras_entrenamiento;

%Armo un vector con todas las muestras a usar en el k-means
Muestras_totales=[Muestras_A_desordenadas([Nro_muestras_entrenamiento+1:Nro_muestras_A],:);Muestras_O_desordenadas([Nro_muestras_entrenamiento+1:Nro_muestras_O],:);Muestras_U_desordenadas([Nro_muestras_entrenamiento+1:Nro_muestras_U],:)];

%Calculo la cantidad de muestras
Nro_total_muestras=length(Muestras_totales);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Grafico de los primeros dos formantes%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Hago el grafico de los primeros 2 formantes para los tres tipos.
%Solo grafico las primeras 40 muestras de entrenamiento.
figure(1)
hold on
plot(Muestras_A_desordenadas([1:Nro_muestras_entrenamiento],1),Muestras_A([1:Nro_muestras_entrenamiento],2),'ro');
plot(Muestras_O_desordenadas([1:Nro_muestras_entrenamiento],1),Muestras_O([1:Nro_muestras_entrenamiento],2),'go');
plot(Muestras_U_desordenadas([1:Nro_muestras_entrenamiento],1),Muestras_U([1:Nro_muestras_entrenamiento],2),'bo');
plot(Media_A(1,1),Media_A(1,2),'rx');
plot(Media_O(1,1),Media_O(1,2),'g*');
plot(Media_U(1,1),Media_U(1,2),'b+');
legend('Formantes a.txt','Formantes o.txt','Formantes u.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%K-MEANS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Armo un vector que tenga la clase de cada muestra
clase_cercana=zeros(1,Nro_total_muestras);
%Condiciones para la primera iteracion
distorsion_paso_previo=0;
distorsion_antigua=0;
distorsion_actualizada=1;
k=2;%Variable que uso para generar figuras.
j=0;%Variable para almacenar la distorsion
%criterio de iteracion:itero hasta que la distorsion varie menos de un 0.001
while(abs((distorsion_actualizada-distorsion_antigua)/distorsion_actualizada)>0.001) && (distorsion_paso_previo~=distorsion_actualizada)
    distorsion_paso_previo=distorsion_antigua;
    distorsion_antigua=distorsion_actualizada;
    distorsion_actualizada=0;
    i=1;%Variable para recorrer todas las muestras en kmeans
    cantidad_clase_1=0;
    cantidad_clase_2=0;
    cantidad_clase_3=0;
    Media_actualizada_A=[0 0];
    Media_actualizada_O=[0 0];
    Media_actualizada_U=[0 0];
    Catalogadas_A=zeros(100,2);
    Catalogadas_O=zeros(100,2);
    Catalogadas_U=zeros(100,2);
%Recorro todas las muestras, clasifico cada muestras en un clase
%dependiendo de la distancia a las medias,suma todos los valores de una
%clase y divido sobre el nro total de muestras de cada clase para obtener
%la media.
   
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
            Catalogadas_A(cantidad_clase_1,:)=Muestras_totales(i,:);
        elseif(clase_cercana(1,i)==2)
            Media_actualizada_O=Media_actualizada_O+Muestras_totales(i,:);
            cantidad_clase_2=cantidad_clase_2+1;
            Catalogadas_O(cantidad_clase_2,:)=Muestras_totales(i,:);
        elseif(clase_cercana(1,i)==3)
            Media_actualizada_U=Media_actualizada_U+Muestras_totales(i,:);
            cantidad_clase_3=cantidad_clase_3+1;
            Catalogadas_U(cantidad_clase_3,:)=Muestras_totales(i,:);
        end
        i=i+1;  
    end

%Dividiendo la suma por la cantidad de cada clase obtengo las nuevas
%medias y divido la distorsion en la cantidad de muestras obtengo la
%distorsion.

Media_actualizada_A=Media_actualizada_A/cantidad_clase_1;
Media_actualizada_O=Media_actualizada_O/cantidad_clase_2;
Media_actualizada_U=Media_actualizada_U/cantidad_clase_3;

j=j+1;
distorsion_actualizada=distorsion_actualizada/Nro_total_muestras;
distorsion_vector(1,j)=distorsion_antigua;
distorsion_vector(2,j)=distorsion_actualizada;

figure(k)
hold on
plot(Catalogadas_A([1:cantidad_clase_1],1),Catalogadas_A([1:cantidad_clase_1],2),'ro');
plot(Catalogadas_O([1:cantidad_clase_2],1),Catalogadas_O([1:cantidad_clase_2],2),'go');
plot(Catalogadas_U([1:cantidad_clase_3],1),Catalogadas_U([1:cantidad_clase_3],2),'bo');
plot(Media_A(1,1),Media_A(1,2),'kx');
plot(Media_O(1,1),Media_O(1,2),'k*');
plot(Media_U(1,1),Media_U(1,2),'k+');
k=k+1;

Media_A=Media_actualizada_A;
Media_O=Media_actualizada_O;
Media_U=Media_actualizada_U;

end

%Grafico la distorsion%
figure(k)
hold on
plot(distorsion_vector(2,:),'r')
legend('Distorsion');

%Calculo los pesos de cada clase
pi_A=cantidad_clase_1/Nro_total_muestras;
pi_O=cantidad_clase_2/Nro_total_muestras;
pi_U=cantidad_clase_3/Nro_total_muestras;

%Calculo las varianzas

Matriz_unos_columnas_A=ones(cantidad_clase_1,1);
Matriz_unos_columnas_O=ones(cantidad_clase_2,1);
Matriz_unos_columnas_U=ones(cantidad_clase_3,1);

Vector_media_A=Matriz_unos_columnas_A*Media_A;
Vector_media_O=Matriz_unos_columnas_O*Media_O;
Vector_media_U=Matriz_unos_columnas_U*Media_U;

Varianza_A=(transpose(Catalogadas_A([1:cantidad_clase_1],:)-Vector_media_A))*(Catalogadas_A([1:cantidad_clase_1],:)-Vector_media_A)/cantidad_clase_1;
Varianza_O=(transpose(Catalogadas_O([1:cantidad_clase_2],:)-Vector_media_O))*(Catalogadas_O([1:cantidad_clase_2],:)-Vector_media_O)/cantidad_clase_2;
Varianza_U=(transpose(Catalogadas_U([1:cantidad_clase_3],:)-Vector_media_U))*(Catalogadas_U([1:cantidad_clase_3],:)-Vector_media_U)/cantidad_clase_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%ACA TERMINO EL ARRANQUE CON KMEANS%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL_antiguo=0;
LL_actualizado=1;
p=1;
while(abs((LL_actualizado-LL_antiguo)/LL_actualizado)>0.001)
    LL_antiguo=LL_actualizado;
    LL_actualizado=0;
    i=1;%Variable para recorrer todas las muestras en kmeans
    cantidad_clase_1=0;
    cantidad_clase_2=0;
    cantidad_clase_3=0;
    Media_actualizada_A=[0 0];
    Media_actualizada_O=[0 0];
    Media_actualizada_U=[0 0];
    Catalogadas_A=zeros(100,2);
    Catalogadas_O=zeros(100,2);
    Catalogadas_U=zeros(100,2);
    phi=zeros(length(Muestras_totales),3);
%Recorro todas las muestras, calculo la responsabilidad, clasifico en base
%a la maxima responsabilidad
    while(i<=Nro_total_muestras)
        %clasifico
        phi(i,:)=funcion_responsabilidad(Muestras_totales(i,:),Media_A,Media_O,Media_U,Varianza_A,Varianza_O,Varianza_U,pi_A,pi_O,pi_U);
        aux = find(max(phi(i,:)) == phi(i,:));
        if(aux==1)
            cantidad_clase_1=cantidad_clase_1+1;
            Catalogadas_A(cantidad_clase_1,:)=Muestras_totales(i,:);
        elseif(aux==2)
            cantidad_clase_2=cantidad_clase_2+1;
            Catalogadas_O(cantidad_clase_2,:)=Muestras_totales(i,:);
        elseif(aux==3)
            cantidad_clase_3=cantidad_clase_3+1;
            Catalogadas_U(cantidad_clase_3,:)=Muestras_totales(i,:);
        end
        i=i+1;
    
    end

    %ACTUALIZO LOS PARAMETROS DE LA MEDIA,VARIANZA Y PESO.
    %Actualizo la media
    j=1;
    phi_clase_1=0;
    phi_clase_2=0;
    phi_clase_3=0;
    while(j<=Nro_total_muestras)
        Media_actualizada_A(1,:)=Media_actualizada_A(1,:)+phi(j,1)*Muestras_totales(j,:);
        phi_clase_1=phi_clase_1+phi(j,1);
        Media_actualizada_O(1,:)=Media_actualizada_O(1,:)+phi(j,2)*Muestras_totales(j,:);
        phi_clase_2=phi_clase_2+phi(j,2);
        Media_actualizada_U(1,:)=Media_actualizada_U(1,:)+phi(j,3)*Muestras_totales(j,:);
        phi_clase_3=phi_clase_3+phi(j,3);
        j=j+1;
    end
    Media_actualizada_A=Media_actualizada_A/phi_clase_1;
    Media_actualizada_O=Media_actualizada_O/phi_clase_2;
    Media_actualizada_U=Media_actualizada_U/phi_clase_3;
    
    %Actualizo la Varianza
    j=1;
    Varianza_actualizada_A=zeros(2,2);
    Varianza_actualizada_O=zeros(2,2);
    Varianza_actualizada_U=zeros(2,2);
    while(j<=Nro_total_muestras)
        Varianza_actualizada_A=Varianza_actualizada_A+phi(j,1)*(transpose(Muestras_totales(j,:)-Media_actualizada_A)*(Muestras_totales(j,:)-Media_actualizada_A));
        Varianza_actualizada_O=Varianza_actualizada_O+phi(j,2)*(transpose(Muestras_totales(j,:)-Media_actualizada_O)*(Muestras_totales(j,:)-Media_actualizada_O));
        Varianza_actualizada_U=Varianza_actualizada_U+phi(j,3)*(transpose(Muestras_totales(j,:)-Media_actualizada_U)*(Muestras_totales(j,:)-Media_actualizada_U));
        j=j+1;
    end
    Varianza_actualizada_A=Varianza_actualizada_A/phi_clase_1;
    Varianza_actualizada_O=Varianza_actualizada_O/phi_clase_2;
    Varianza_actualizada_U=Varianza_actualizada_U/phi_clase_3;
    %Actualizo los pi_k
    pi_A_actualizado=phi_clase_1/Nro_total_muestras;
    pi_O_actualizado=phi_clase_2/Nro_total_muestras;
    pi_U_actualizado=phi_clase_3/Nro_total_muestras;
    
    %Actualizo LL%%%%%%%
    j=1;
    LL_actualizado=0;
    while(j<=Nro_total_muestras)
        p1=mvnpdf(Muestras_totales(j,:),Media_actualizada_A,Varianza_actualizada_A)*pi_A_actualizado;
        p2=mvnpdf(Muestras_totales(j,:),Media_actualizada_O,Varianza_actualizada_O)*pi_O_actualizado;
        p3=mvnpdf(Muestras_totales(j,:),Media_actualizada_U,Varianza_actualizada_U)*pi_U_actualizado;
        LL_actualizado=LL_actualizado+log(p1+p2+p3);
        j=j+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Actualizo las variables para la prox
    %%%%%%%%%%%%%%%%%%%%%%%%%%iteracion%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Media_A=Media_actualizada_A;
    Media_O=Media_actualizada_O;
    Media_U=Media_actualizada_U;
    Varianza_A=Varianza_actualizada_A;
    Varianza_O=Varianza_actualizada_O;
    Varianza_U=Varianza_actualizada_U;
    pi_A=pi_A_actualizado;
    pi_O=pi_O_actualizado;
    pi_U=pi_U_actualizado;
    
    
   p=p+1; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%PRUEBA%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Muestras_evaluacion=[Muestras_A_desordenadas([Nro_muestras_A+1:end],:);Muestras_O_desordenadas([Nro_muestras_O+1:end],:);Muestras_U_desordenadas([Nro_muestras_U+1:end],:)];
Longitud_test=length(Muestras_evaluacion);
i=1;
while(i<=Longitud_test)
    %Agrego un 1 en log pk
    GKA=(-0.5)*log(det(Varianza_A))-0.5*(Muestras_evaluacion(i,:)-Media_A)*(Varianza_A^-1)*transpose(Muestras_evaluacion(i,:)-Media_A)+log(pi_A);
    GKO=(-0.5)*log(det(Varianza_O))-0.5*(Muestras_evaluacion(i,:)-Media_O)*(Varianza_O^-1)*transpose(Muestras_evaluacion(i,:)-Media_O)+log(pi_O);
    GKU=(-0.5)*log(det(Varianza_U))-0.5*(Muestras_evaluacion(i,:)-Media_U)*(Varianza_U^-1)*transpose(Muestras_evaluacion(i,:)-Media_U)+log(pi_U);

   valores=[GKA GKO GKU];
   respuesta_2(i) = find(max(valores) == valores);
   i=i+1;
end

figure(k+1)
hold on
plot(Muestras_evaluacion(:,1),Muestras_evaluacion(:,2),'ko')
plot(Media_A(1,1),Media_A(1,2),'rx');
plot(Media_O(1,1),Media_O(1,2),'g*');
plot(Media_U(1,1),Media_U(1,2),'b+');
