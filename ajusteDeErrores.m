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

function [solucion] = aproximacionExponencial(listX,listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1,5);
  for i=1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = listX(1,i).**2;
    matrizAproximacion(i,4) = reallog(listY(1,i));
    matrizAproximacion(i,5) = listX(1,i)*reallog(listY(1,i));
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
  
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,5); matrizAproximacion(cantFilas(1,2)+1,4)] ;
  
  Solucion = Matrix1\Matrix2;
  Solucion(2,1) = exp(Solucion(2,1));
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
   ok=0;
   while(ok == 0)
     seleccionMenuOpciones = inputdlg({"Elija la accion que desee realizar:\n\n1- Aproximar\n\n2- Comparar aproximaciones\n\n3- Finalizar\n\n"},"Ajuste App", [0.5]);
     cargoOk=0
     switch(str2num(seleccionMenuOpciones{1}))
       case(1)%Aproximacion Lineal
         while(cargoOk == 0)
           cargoOk = interfazAjuste()
         endwhile
       case(2)
           interfazComparaciones() 
       case(3)
           ok = 1;
           exit;
       otherwise
            errordlg("La opcion requerida es incorrecta, por favor vuelvalo a intentar", "Error de peticion");
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

function [x] = valores(parametro)
  x = inputdlg({cstrcat("Ingrese una lista de valores para ",parametro)},"Ajuste App",[0.5]);
  while (isempty(str2num(x{1}))) 
      errordlg(cstrcat("No ingreso los valores de ",parametro), "Error al procesar");
      x = inputdlg({cstrcat("Ingrese una lista de valores para ",parametro)},"Ajuste App",[0.5]);
  endwhile
endfunction

function esOK = interfazAjuste()
    datosAproximacion = inputdlg({"Ingrese la lista de valores para X","Ingrese la lista de valores de Y", "Ingrese la cantidad de decimales de la aproximacion"},"Ajuste App",[0.5]);
    if(or(isempty(datosAproximacion{1}), isempty(datosAproximacion{2}), isempty(datosAproximacion{3})))
      errordlg("Ha ingresado incorrectamente los datos para generar la aproximacion", "Error al procesar");
      esOK = 0;
    else
      listX = str2num(datosAproximacion{1});
      listY = str2num(datosAproximacion{2});
      decimales = str2num(datosAproximacion{3});
      esOK = 1;
      eligioBien = 0;
      while(eligioBien == 0)
      %Me falta manejar el error
      seleccionMenuOpcionesMostrar = inputdlg({"Elija opcion de ajuste:\n\n-1-Recta de minimos cuadrados\n\n-2-Parabola de minimos cuadrados\n\n-3-Aproximacion Exponencial\n\n-4- Aproximacion Potencial\n\n-5-Aproximacion Hiperbola\n\n"},"Ajuste App", [0.5]);
          switch(str2num(seleccionMenuOpcionesMostrar{1}))
               case(1)%Aproximacion Lineal 
                 [Solucion]= aproximacionLineal(listX,listY);
                 msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(1,1))," x+ ", num2str(Solucion(2,1))),"Ajuste App");
                 eligioBien = 1;
              case(2)
                   [Solucion]= aproximacionParabola(listX,listY);
                   msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(1,1))," x^2+ ",num2str(Solucion(2,1))," x ",num2str(Solucion(3,1))),"Ajuste App");
                   eligioBien = 1;
              case(3)                   
                 [Solucion]=aproximacionExponencial(listX,listY);
                 msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(2,1))," e ^ ",num2str(Solucion(1,1)))," x ","Ajuste App");
                  eligioBien = 1;
              case(4)
                  [Solucion] = aproximacionPotencial(listX,listY);
                  msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(2,1))," x ^ ",num2str(Solucion(1,1))),"Ajuste App");
                  eligioBien = 1;
              case(5)
                  [Solucion]=aproximacionHiperbolica(listX,listY);
                  msgbox(cstrcat("La funcion es:\n","y = ",num2str(Solucion(2,1))," x ^ ",num2str(Solucion(1,1))),"Ajuste App");
                  eligioBien = 1;
              otherwise
                 errordlg("Error al elegir el tipo de aproximacion, por favor vuelca a intentarlo"
            endswitch
        endwhile
     endif 
endfunction
















