%--Carga de packages



function main

  interfazPrincipal()
endfunction
%Auxiliares
function [Solucion] = aproximacionPotencial(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 6);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i, 3) = reallog(listX(1,i));
    matrizAproximacion(i, 4) = reallog(listX(1,i)).**2;
    matrizAproximacion(i,5) = reallog(listY(1,i));
    matrizAproximacion(i,6) = reallog(listX(1,i))*reallog(listY(1,i));
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6)) ; 
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,3) ; matrizAproximacion(cantFilas(1,2)+1,3), cantFilas(1,2)]
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,6); matrizAproximacion(cantFilas(1,2)+1,5)];
  Solucion = Matrix1\Matrix2;
  sizeMatrix2 = size(Matrix2);
  Solucion(sizeMatrix2(1,1),1) = exp(Solucion(sizeMatrix2(1,1),1));
endfunction

function [Solucion] = aproximacionLineal(listX, listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = listX(1,i).**2;
    matrizAproximacion(i,4) = listX(1,i)*listY(1,i);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,4); matrizAproximacion(cantFilas(1,2)+1,2)] ;
  Solucion = Matrix1\Matrix2;
 
endfunction

function [Solucion] = aproximacionHiperbolica(listX, listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i, 3) = listX(1,i).**2;
    matrizAproximacion(i, 4) = listX(1,i)*listY(1,i);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  Matrix1 = [cantFilas(1,2),matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1),matrizAproximacion(cantFilas(1,2)+1,3)];
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,2);matrizAproximacion(cantFilas(1,2)+1,4)];
  Solucion = Matrix1\Matrix2;
  Solucion(1,1) = Solucion(1,1)*(Solucion(2,1).**(-1));
  Solucion(2,1) = Solucion(2,1).**(-1);
  
 
endfunction

function [Solucion] = aproximacionParabola(listX, listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 7);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = listX(1,i).**2;
    matrizAproximacion(i,4) = listX(1,i).**3;
    matrizAproximacion(i,5) = listX(1,i).**4;
    matrizAproximacion(i,6) = listX(1,i)*listY(1,i);
    matrizAproximacion(i,7) = listY(1,i)*listX(1,i).**2;
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6));
  matrizAproximacion(cantFilas(1,2)+1,7) = sum(matrizAproximacion(:,7));
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,5);matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4);cantFilas(1,2), matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3)];
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,7); matrizAproximacion(cantFilas(1,2)+1,6);matrizAproximacion(cantFilas(1,2)+1,2)];
  Solucion = Matrix1\Matrix2;
  
endfunction

function interfazPrincipal()
   pkg load control
   pkg load symbolic
   seleccionMenuOpciones = inputdlg({"Elija la accion que desee realizar:\n\n-1-Aproximar\n\n-2-Comparar aproximaciones\n\n -3-Finalizar\n\n"},"Ajuste App", [0.5]);
   ok=0;
while(ok == 0)
      switch(str2num(seleccionMenuOpciones{1}))
        case(1)%Aproximacion Lineal
          interfazAjuste()
        case(2)
          interfazComparaciones() 
        case(3)
           ok = 1;
           exit;
         endswitch
  endwhile
endfunction

function ploteoPuntos(listX,listY)
  x = listX;
  y = listY;
  plot(x,y,'rx','MarkerSize',8); 
  xlabel('x');
  ylabel('y');
endfunction

function interfazComparaciones()
  endfunction

function [x] valores(parametro)
  x = inputdlg({cstrcat("Ingrese una lista de valores para ",parametro)},"Ajuste App",[0.5]);
  while (isempty(str2num(auxiliarX{1}))) 
      errordlg(cstrcat("No ingreso los valores de ",parametro), "Error al procesar");
      x = inputdlg({cstrcat("Ingrese una lista de valores para ",parametro)},"Ajuste App",[0.5]);
  endwhile
  

function interfazAjuste()
  [auxiliarX] = valores("x");
  [auxiliarY] = valores("y")

  listX = str2num(auxiliarX{1});
  listY = str2num(auxiliarY{1});

    seleccionMenuOpcionesMostrar = inputdlg({"Elija opcion de ajuste:\n\n-1-Recta de minimos cuadrados\n\n-2-Parabola de minimos cuadrados\n\n-3-Aproximacion Exponencial\n\n-4- Aproximacion Potencial\n\n-5-Aproximacion Hiperbola\n\n"},"Ajuste App", [0.5]);
    ok=0;
       switch(str2num(seleccionMenuOpcionesMostrar{1}))
            case(1)%Aproximacion Lineal 
              [Solucion]= aproximacionLineal(listX,listY);
              msgbox(cstrcat("Se utilizó el metodo: Recta de minimos cuadrados \n\nLa funcion es:\n","y = ",num2str(Solucion(1,1))," x+ ", num2str(Solucion(2,1))),"Ajuste App");
              ploteoPuntos(listX,listY);
           case(2)
              [Solucion]= aproximacionParabola(listX,listY);
              msgbox(cstrcat("Se utilizó el metodo: Parabola de minimos cuadrados \n\nLa funcion es:\n","y = ",num2str(Solucion(1,1))," x^2+ ",num2str(Solucion(2,1))," x ",num2str(Solucion(3,1))),"Ajuste App");
               ploteoPuntos(listX,listY);
           case(3)
               msgbox("Se utilizó el metodo: Aproximacion Exponencial\n\n");
              %aproximacionExponencial(listX,listY);
           case(4)
             [Solucion] = aproximacionPotencial(listX,listY);
             msgbox(cstrcat("Se utilizó el metodo: Aproximacion Potencial\n\nLa funcion es:\n","y = ",num2str(Solucion(2,1))," x ^ ",num2str(Solucion(1,1))),"Ajuste App");
           case(5)
            [Solucion]= aproximacionHiperbolica(listX,listY);
              msgbox("Se utilizó el metodo: Aproximacion Hiperbola \n\n");
             %msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(2,1))," x ^ ",num2str(Solucion(1,1))),"Ajuste App");
         endswitch
   uiwait();
  interfazPrincipal();
endfunction
















