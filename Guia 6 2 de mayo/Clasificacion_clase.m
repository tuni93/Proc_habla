function [clase,distancia] = Clasificacion_clase(x,media_1,media_2,media_3)
        D1=sqrt((x(1,1)-media_1(1,1))^2+(x(1,2)-media_1(1,2))^2);
        D2=sqrt((x(1,1)-media_2(1,1))^2+(x(1,2)-media_2(1,2))^2);
        D3=sqrt((x(1,1)-media_3(1,1))^2+(x(1,2)-media_3(1,2))^2);

        if(D1<=D2) && (D1<=D3)
            clase=1;
            distancia=D1;
        elseif (D2<=D1) && (D2<=D3)
            clase=2;
            distancia=D2;
        elseif (D3<=D1) && (D3<=D2)     
            clase=3;
            distancia=D3;
                     
end
