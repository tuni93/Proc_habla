%Guia 6
%Funes Pablo Nicolas
%Padron 94894


%Lectura de los archivos

Muestras_A=load('a.txt');
Muestras_O=load('o.txt');
Muestras_U=load('u.txt');

%Tomo parte de las muestras para entrenamiento.
Nro_muestras_entrenamiento=40;
%Elijo los primeros Nro_formantes.
Nro_formantes=2;

%Selecciono en base a las muestras y formantes, las muestras de
%entrenamiento

Muestras_entrenamiento_A=Muestras_A([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
Muestras_entrenamiento_O=Muestras_O([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
Muestras_entrenamiento_U=Muestras_U([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);

%Hago las elipses
puntos_elipse=0:1:49;
puntos_elipse=puntos_elipse*2*pi/49;
elipse=zeros(2,50);
elipse(1,:)=sin(puntos_elipse);
elipse(2,:)=cos(puntos_elipse);



%Hago el grafico de los primeros 2 formantes para los tres tipos.
figure(1)
hold on
plot(Muestras_entrenamiento_A(:,1),Muestras_entrenamiento_A(:,2),'ro');
plot(Muestras_entrenamiento_O(:,1),Muestras_entrenamiento_O(:,2),'go');
plot(Muestras_entrenamiento_U(:,1),Muestras_entrenamiento_U(:,2),'bo');
legend('Formantes a.txt','Formantes o.txt','Formantes u.txt');


%Calculo la media de las muestras de entrenamiento para cada tipo.
Matriz_unos=ones(1,Nro_muestras_entrenamiento);
%Suma de todos los elementos.
Media_A=Matriz_unos*Muestras_entrenamiento_A;
Media_O=Matriz_unos*Muestras_entrenamiento_O;
Media_U=Matriz_unos*Muestras_entrenamiento_U;


%Divido por el nro de muestras obtengo el promedio.
Media_A=Media_A/Nro_muestras_entrenamiento;
Media_O=Media_O/Nro_muestras_entrenamiento;
Media_U=Media_U/Nro_muestras_entrenamiento;


%Calculo la varianza.
Matriz_unos_columnas=ones(Nro_muestras_entrenamiento,1);
Vector_media_A=Matriz_unos_columnas*Media_A;
Vector_media_O=Matriz_unos_columnas*Media_O;
Vector_media_U=Matriz_unos_columnas*Media_U;
%Aca hago hago la multiplicacion de muestras menos la media;
Varianza_A=(transpose(Muestras_entrenamiento_A-Vector_media_A))*(Muestras_entrenamiento_A-Vector_media_A);
Varianza_O=(transpose(Muestras_entrenamiento_O-Vector_media_O))*(Muestras_entrenamiento_O-Vector_media_O);
Varianza_U=(transpose(Muestras_entrenamiento_U-Vector_media_U))*(Muestras_entrenamiento_U-Vector_media_U);
%Me tengo que quedar solo con los elementos varianza(1,1) y varianza(2,2).
Var_A=[Varianza_A(1,1),Varianza_A(2,2)];
Var_O=[Varianza_O(1,1),Varianza_O(2,2)];
Var_U=[Varianza_U(1,1),Varianza_U(2,2)];
%Divido por el nro de muestras
Var_A=Var_A/Nro_muestras_entrenamiento;
Var_O=Var_O/Nro_muestras_entrenamiento;
Var_U=Var_U/Nro_muestras_entrenamiento;

hold on
plot(Media_A(1,1)+sqrt(Var_A(1,1))*elipse(1,:),Media_A(1,2)+sqrt(Var_A(1,2))*elipse(2,:),'r');
plot(Media_O(1,1)+sqrt(Var_O(1,1))*elipse(1,:),Media_O(1,2)+sqrt(Var_O(1,2))*elipse(2,:),'g');
plot(Media_U(1,1)+sqrt(Var_U(1,1))*elipse(1,:),Media_U(1,2)+sqrt(Var_U(1,2))*elipse(2,:),'b');





%Calculo el promedio de las varianzas
Varianza=Varianza_A+Varianza_O+Varianza_U;
Varianza=Varianza/3;


%Para verifica a que categoria o clase pertenece uso la funcion g
%Elijo la muestras para analizar
%Muestra_analisis=transpose(Muestras_U(45,[1:Nro_formantes]));

%Armo un vector de muestras para probar el rendimiento.
Muestras_prueba=[Muestras_A([41:50],[1:Nro_formantes]);Muestras_O([41:50],[1:Nro_formantes]);Muestras_U([41:50],[1:Nro_formantes])];
Muestras_prueba_t=transpose(Muestras_prueba);

i=1;
while(i<=30)

Muestra_analisis=Muestras_prueba_t(:,i);    
XOk_A=Media_A*((Varianza)^-1)*transpose(Media_A);
WOk_A=transpose(((Varianza)^-1)*transpose(Media_A));
GKA=(-1/2)*XOk_A+WOk_A*Muestra_analisis;

XOk_O=Media_O*((Varianza)^-1)*transpose(Media_O);
WOk_O=transpose(((Varianza)^-1)*transpose(Media_O));
GKO=(-1/2)*XOk_O+WOk_O*Muestra_analisis;

XOk_U=Media_U*((Varianza)^-1)*transpose(Media_U);
WOk_U=transpose(((Varianza)^-1)*transpose(Media_U));
GKU=(-1/2)*XOk_U+WOk_U*Muestra_analisis;

%1 clase A 2 clase O 3 clase U

if(GKA>GKO)&&(GKA>GKU)
    respuesta=1;
end
if(GKO>GKA)&&(GKO>GKA)
    respuesta=2;
end
if(GKU>GKO) && (GKU>GKA)
    respuesta=3;
end
Resultados(i)=respuesta;
i=i+1;
end



