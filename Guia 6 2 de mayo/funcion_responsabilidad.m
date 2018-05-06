function [phi] = funcion_responsabilidad(X,Media_1,Media_2,Media_3,Varianza_1,Varianza_2,Varianza_3,pi_1,pi_2,pi_3)
 
     %Calculo la probabilidad total
     p1=mvnpdf(X,Media_1,Varianza_1);
     p2=mvnpdf(X,Media_2,Varianza_2);
     p3=mvnpdf(X,Media_3,Varianza_3);
     
     Ptotal=p1*pi_1+p2*pi_2+p3*pi_3;
   
     phi(1,1)=p1*pi_1/Ptotal;
     phi(1,2)=p2*pi_2/Ptotal;
     phi(1,3)=p3*pi_3/Ptotal;
end
