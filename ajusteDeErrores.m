%--Carga de packages
function main
  interfazPrincipal()
endfunction

%Auxiliares
function [funcion] = aproximacionPotencial(listX, listY)
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
  
  funcion = @(x)Solucion(2,1)*x.^(Solucion(1,1));
  
endfunction

function [funcion] = aproximacionExponencial(listX,listY)
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
  
   funcion = @(x)Solucion(2,1)*e.^(Solucion(1,1)*x) ;
endfunction 

function [funcion] = aproximacionLineal(listX, listY)
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
  
  funcion = @(x)x*Solucion(1,1)+Solucion(2,1);
 
endfunction

function [funcion] = aproximacionHiperbolica(listX, listY)
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
    
    funcion = @(x)Solucion(1,1)*((x+Solucion(2,1)).^(-1));
   
endfunction

function [funcion] = aproximacionParabola(listX, listY)
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
  
  funcion = @(x)x.^2*Solucion(1,1)+x*Solucion(2,1)+Solucion(3,1);
  
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
           return;
       otherwise
            errordlg("La opcion requerida es incorrecta, por favor vuelvalo a intentar", "Error de peticion");
    endswitch
  endwhile
endfunction

function ploteoPuntos(listX,listY,funcion)
  x = listX;
  y = listY;
  hold on;
  plot(x,y,'rx','MarkerSize',8); 
  a = min(listX);
  b = max(listX);
  fplot(funcion,[a b]);
  xlabel('x');
  ylabel('y');
  legend('Puntos', 'Curva representativa');
  uiwait;
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

function [seleccion, flag] = mostrarMenuSeleccionFuncionalidad()
  listaDeOpciones = {"Mostrar funcion aproximante", "Detalles de calculos", "Grafica"};
  [seleccion, flag] = listdlg("ListString", listaDeOpciones, "PromptString", "Seleccione la funcionalidad deseada:", "Name", "Menu Funcionalidad", "ListSize", [290 140], "SelectionMode", "Single", "OKString", "Ok", "CancelString", "Cancelar");
endfunction

function [tipoFuncion, flag] = mostrarMenuFuncionesAproximantes()
  listaDeOpciones = {"Recta de minimos cuadrados", "Parabola de minimos cuadrados", "Aproximacion Potencial", "Aproximacion Exponencial", "Aproximacion Hiperbolica"};
  [tipoFuncion, flag] = listdlg("ListString", listaDeOpciones, "PromptString", "Seleccione la funcion que desea utilizar para aproximar:", "Name", "Menu Funcion", "ListSize", [290 140], "SelectionMode", "Single", "OKString", "Ok", "CancelString", "Cancelar");
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
      [tipoFuncion, flagFuncion] = mostrarMenuFuncionesAproximantes();
      if(flagFuncion == 1)
        bandera = 1;
        while(bandera == 1)
          [caracteristica, bandera] = mostrarMenuSeleccionFuncionalidad();
           switch(caracteristica)
             case(1)
               mostrarFuncionAproximante(tipoFuncion, listX, listY);
             case(2)
               mostrarDetalleDeCalculo(tipoFuncion, listX, listY);
             case(3)
               graficar(tipoFuncion, listX, listY);
           endswitch
          endwhile
      endif
     endif
endfunction

function mostrarDetalleDeCalculo(tipoFuncion, listX, listY)
endfunction

function graficar(tipoFuncion, listX, listY)
  [funcion]= obtenerFuncion(tipoFuncion,listX,listY);
  ploteoPuntos(listX,listY,funcion);
endfunction

function mostrarFuncionAproximante(tipoFuncion, listX, listY)
  [funcion]= obtenerFuncion(tipoFuncion,listX,listY);
  syms x
  g = symfun(funcion,x);
  msgbox(cstrcat("La funcion es: \n\n",disp(g),"\n\n"));
endfunction

function [funcion] = obtenerFuncion(tipoFuncion,listX,listY)
  switch(tipoFuncion)
    case(1)
      [funcion]= aproximacionLineal(listX,listY);
    case(2)
      [funcion]= aproximacionParabola(listX,listY);
    case(3)
       [funcion]=aproximacionExponencial(listX,listY);
    case(4)
      [funcion] = aproximacionPotencial(listX,listY);
    case(5)
      [funcion]=aproximacionHiperbolica(listX,listY);
  endswitch
endfunction








