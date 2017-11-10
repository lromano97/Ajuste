%--Carga de packages
function main
  interfazPrincipal()
endfunction

function valor = redondear(numero, decimales)
  if(numero > -1 && numero < 1)
    valor = chop(numero, decimales, 10);
  else
    valor = chop(numero, decimales+1, 10);
  endif
endfunction

%Auxiliares
function [matrizAproximacion] = aproximacionPotencial(listX, listY, decimales)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 6);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i, 3) = redondear(reallog(listX(1,i)), decimales);
    matrizAproximacion(i, 4) = redondear(reallog(listX(1,i)).**2, decimales);
    matrizAproximacion(i,5) = redondear(reallog(listY(1,i)), decimales);
    matrizAproximacion(i,6) = redondear(reallog(listX(1,i))*reallog(listY(1,i)), decimales);
  endfor
  mat2str(matrizAproximacion, decimales+1);
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6)); 
endfunction

function [matrizAproximacion] = aproximacionExponencial(listX,listY, decimales)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1,5);
  for i=1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = redondear(listX(1,i).**2, decimales);
    matrizAproximacion(i,4) = redondear(reallog(listY(1,i)), decimales);
    matrizAproximacion(i,5) = redondear(listX(1,i)*reallog(listY(1,i)), decimales);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
endfunction 

function [matrizAproximacion] = aproximacionLineal(listX, listY, decimales)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = redondear(listX(1,i).**2, decimales);
    matrizAproximacion(i,4) = redondear(listX(1,i)*listY(1,i), decimales);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
endfunction

function [matrizAproximacion] = aproximacionHiperbolica(listX, listY, decimales)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i, 3) = redondear(listX(1,i).**2, decimales);
    matrizAproximacion(i, 4) = redondear(listX(1,i)*listY(1,i), decimales);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
endfunction

function [matrizAproximacion] = aproximacionParabola(listX, listY, decimales)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,2)+1, 7);
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i);
    matrizAproximacion(i,2) = listY(1,i);
    matrizAproximacion(i,3) = redondear(listX(1,i).**2, decimales);
    matrizAproximacion(i,4) = redondear(listX(1,i).**3, decimales);
    matrizAproximacion(i,5) = redondear(listX(1,i).**4, decimales);
    matrizAproximacion(i,6) = redondear(listX(1,i)*listY(1,i), decimales);
    matrizAproximacion(i,7) = redondear(listY(1,i)*listX(1,i).**2, decimales);
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1));
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2));
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3));
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4));
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5));
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6));
  matrizAproximacion(cantFilas(1,2)+1,7) = sum(matrizAproximacion(:,7));
endfunction

function interfazPrincipal()
   pkg load control
   pkg load symbolic
   ok=0;
   while(ok == 0)
    Funcionalidades = {"Aproximar", "Comparar aproximaciones", "Finalizar"};
    [Funcionalidad, flag] = listdlg("ListString", Funcionalidades, "PromptString", "Seleccione la funcionalidad:", "Name", "Menu Funcion", "ListSize", [290 140], "SelectionMode", "Single", "OKString", "Ok", "CancelString", "Cancelar");
      if(flag == 1)
       cargoOk=0
       switch(Funcionalidad)
         case(1)%Aproximacion Lineal
           while(cargoOk == 0)
             cargoOk = interfazAjuste()
           endwhile
         case(2)
             interfazComparaciones() 
         case(3)
             respuesta = questdlg("Esta seguro que desea salir de AMIC?", "AMIC", "Si", "No", "No");
             if(strcmp(respuesta, "Si"))
              msgbox("Gracias por utilizar AMIC", "AMIC", "none");
              ok = 1;
              return;
           endif
        endswitch
       else
        respuesta = questdlg("Esta seguro que desea salir de AMIC?", "AMIC", "Si", "No", "No");
        if(strcmp(respuesta, "Si"))
          msgbox("Gracias por utilizar AMIC", "AMIC", "none");
          ok = 1;
          return;
        endif
       endif
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

function [Error] = cargarValoresXeY(cantidadElecciones,listX,listY)
  cantFilas = size(listX);
  Error = zeros(cantFilas(1,2),3+cantidadElecciones+cantidadElecciones);
  for i = 1:cantFilas(1,2)
    Error(i,1) = i;
    Error(i,2) = listX(1,i);
    Error(i,3) = listY(1,i);
    endfor
  endfunction
  
function errorMinimo(compararErrores, decimales)
  cantidadAproxs = size(compararErrores);
  j = 1;
  for i = 1:cantidadAproxs(1,2)
    if(compararErrores(1,i)!=0)
    matrizMinima(1,j) = compararErrores(1,i);
    j = j+1;
    endif
    endfor
  valorMinimo = min(matrizMinima);
  switch(valorMinimo)
    case(compararErrores(1,1))
      msgbox(cstrcat("La recta con: ", mat2str(redondear(valorMinimo, decimales))), "Menor error:\n\n" , "none");
    case(compararErrores(1,2))
      msgbox(cstrcat("La parabola con: ", mat2str(redondear(valorMinimo, decimales))), "Menor error:\n\n" , "none");
    case(compararErrores(1,3))
      msgbox(cstrcat("La potencial con: ", mat2str(redondear(valorMinimo, decimales))), "Menor error:\n\n" , "none");
    case(compararErrores(1,4))
      msgbox(cstrcat("La exponencial con: ", mat2str(redondear(valorMinimo, decimales))), "Menor error:\n\n" , "none");
    case(compararErrores(1,5))
      msgbox(cstrcat("La hiperbolica con: ", mat2str(redondear(valorMinimo, decimales))), "Menor error:\n\n" , "none");
  endswitch
endfunction

function interfazComparaciones()
  datosAproximacion = inputdlg({"Ingrese la lista de valores para X","Ingrese la lista de valores de Y", "Ingrese la cantidad de decimales de la aproximacion"},"Ajuste App",[0.5]);
  if(or(isempty(datosAproximacion{1}), isempty(datosAproximacion{2}), isempty(datosAproximacion{3})))
      errordlg("Ha ingresado incorrectamente los datos para generar la aproximacion", "Error al procesar");
      esOK = 0;
    else
      listX = str2num(datosAproximacion{1});
      listY = str2num(datosAproximacion{2});
      decimales = str2num(datosAproximacion{3});
      cantFilas = size(listX);
      listaDeOpciones = {"Recta de minimos cuadrados", "Parabola de minimos cuadrados", "Aproximacion Potencial", "Aproximacion Exponencial", "Aproximacion Hiperbolica"};
      [Eleccion,ok] = listdlg ("ListString", listaDeOpciones,"SelectionMode", "Multiple"); 
      cantidadElecciones = numel(Eleccion);
      [Error] = cargarValoresXeY(cantidadElecciones,listX,listY);
      compararErrores = zeros(1,5);
      if (ok == 1)
        for i = 1:cantidadElecciones
        switch(Eleccion(i))
        case(1)
          [matrizAproximacion] = aproximacionLineal(listX, listY, decimales);
          Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
          Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,4); matrizAproximacion(cantFilas(1,2)+1,2)] ;
          Solucion = Matrix1\Matrix2;
          Solucion(1,1) = redondear(Solucion(1,1),decimales);
          Solucion(2,1) = redondear(Solucion(2,1),decimales);
          a = Solucion(1,1);
          b = Solucion(2,1);
          for j = 1:cantFilas(1,2)
          Error(j,i+3+cantidadElecciones) = ((listY(:,j)-(listX(:,j)*a+b)).**2);
          compararErrores(1,1) += ((listY(:,j)-(listX(:,j)*a+b)).**2);
          endfor
        case(2)
          [matrizAproximacion] = aproximacionParabola(listX, listY, decimales);
          Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,5);matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4);cantFilas(1,2), matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3)];
          Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,7); matrizAproximacion(cantFilas(1,2)+1,6);matrizAproximacion(cantFilas(1,2)+1,2)];
          Solucion = Matrix1\Matrix2;
          Solucion(1,1) = redondear(Solucion(1,1),decimales);
          Solucion(2,1) = redondear(Solucion(2,1),decimales);
          a = Solucion(1,1);
          b = Solucion(2,1);
          c = Solucion(3,1);
          for j = 1:cantFilas(1,2)
           Error(j,i+3+cantidadElecciones) = (listY(:,j)-(((listX(:,j).**2)*a)+((listX(:,j))*b)+c)).**2;
           compararErrores(1,2) += (listY(:,j)-(((listX(:,j).**2)*a)+((listX(:,j))*b)+c)).**2;
          endfor
        case(3)
          [matrizAproximacion] = aproximacionPotencial(listX, listY, decimales);
          Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,3) ; matrizAproximacion(cantFilas(1,2)+1,3), cantFilas(1,2)]
          Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,6); matrizAproximacion(cantFilas(1,2)+1,5)];
          Solucion = Matrix1\Matrix2;
          sizeMatrix2 = size(Matrix2);
          Solucion(1,1) = redondear(Solucion(1,1),decimales);
          Solucion(sizeMatrix2(1,1),1) = redondear(exp(Solucion(sizeMatrix2(1,1),1)),decimales);
          Solucion(2,1) = redondear(Solucion(2,1),decimales);
          a = Solucion(2,1);
          b = Solucion(1,1);
          for j = 1:cantFilas(1,2)
            Error(j,i+3+cantidadElecciones) = (listY(:,j)-(a*(listX(:,j).**b))).**2;
            compararErrores(1,3) += (listY(:,j)-(a*(listX(:,j).**b))).**2;
          endfor
         case(4)
          [matrizAproximacion] = aproximacionExponencial(listX, listY, decimales);
          Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
          Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,5); matrizAproximacion(cantFilas(1,2)+1,4)] ;
          Solucion = Matrix1\Matrix2;
          Solucion(1,1) = redondear(Solucion(1,1),decimales);
          Solucion(2,1) = redondear(exp(Solucion(2,1)),decimales);
          a = Solucion(1,1);
          b = Solucion(2,1);
          for j = 1:cantFilas(1,2)
            Error(j,i+3+cantidadElecciones) = (listY(:,j)-(b*(e.**(a*listX(:,j))))).**2; 
            compararErrores(1,4) += (listY(:,j)-(b*(e.**(a*listX(:,j))))).**2; 
          endfor
          case(5)
          [matrizAproximacion] = aproximacionHiperbolica(listX, listY, decimales);
          Matrix1 = [cantFilas(1,2),matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1),matrizAproximacion(cantFilas(1,2)+1,3)];
          Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,2);matrizAproximacion(cantFilas(1,2)+1,4)];
          Solucion = Matrix1\Matrix2;
          Solucion(1,1) = redondear(Solucion(1,1)*(Solucion(2,1).**(-1)),decimales);
          Solucion(2,1) = redondear(Solucion(2,1).**(-1),decimales);
          a = Solucion(1,1);
          b = Solucion(2,1);
          for j = 1:cantFilas(1,2)
            Error(j,i+3+cantidadElecciones) = (listY(:,j)-(b*((listX(:,j)+a).**(-1)))).**2;
            compararErrores(1,5) += (listY(:,j)-(b*((listX(:,j)+a).**(-1)))).**2;
          endfor
          endswitch
        endfor
        msgbox(disp(Error), "Comparacion de aproximaciones:\n\n", "none");
       errorMinimo(compararErrores, decimales); 
      endif
   endif
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
               mostrarFuncionAproximante(tipoFuncion, listX, listY, decimales);
             case(2)
               mostrarDetalleDeCalculo(tipoFuncion, listX, listY, decimales);
             case(3)
               graficar(tipoFuncion, listX, listY, decimales);
           endswitch
          endwhile
      endif
     endif
endfunction

function mostrarDetalleDeCalculo(tipoFuncion, listX, listY, cantidadDecimales)
  cantFilas = size(listX);
  switch(tipoFuncion)
    case(1)
      [matrizAproximacion] = aproximacionLineal(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,4); matrizAproximacion(cantFilas(1,2)+1,2)] ;
      msgbox(disp(matrizAproximacion), "Detalle de calculos:\n\n", "none"); 
      msgbox(cstrcat(mat2str(Matrix1(1,1)), "*a", "+ b*", mat2str(Matrix1(1,2)), " = ", mat2str(Matrix2(1,1)), "\n\n", mat2str(Matrix1(2,1)), "* a + ", mat2str(Matrix1(2,2)), " * b = ", mat2str(Matrix2(2,1))), "Sistema utilizado", "none");
    case(2)
      [matrizAproximacion] = aproximacionParabola(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,5);matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4);cantFilas(1,2), matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,7); matrizAproximacion(cantFilas(1,2)+1,6);matrizAproximacion(cantFilas(1,2)+1,2)];
      msgbox(disp(matrizAproximacion), "Detalle de calculos:\n\n", "none"); 
      msgbox(cstrcat(mat2str(Matrix1(1,1)), " * a + ", mat2str(Matrix1(1,2)), " * B = ", mat2str(Matrix2(1,1)), "\n\n", mat2str(Matrix1(2,1)), " * a + ", mat2str(Matrix1(2,2)), " * B = ", mat2str(Matrix2(2,1))), "Sistema utilizado" , "none"); 
    case(3)
      [matrizAproximacion] = aproximacionPotencial(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,3) ; matrizAproximacion(cantFilas(1,2)+1,3), cantFilas(1,2)]
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,6); matrizAproximacion(cantFilas(1,2)+1,5)];
      msgbox(disp(matrizAproximacion), "Detalle de calculos:\n\n", "none"); 
      msbox(cstrcat(mat2str(Matrix1(1,1)), " * a + ", mat2str(Matrix1(1,2)), " * B = ", mat2str(Matrix2(1,1)), "\n\n", mat2str(Matrix1(2,1)), " * a + ", mat2str(Matrix1(2,2)), " * B = ", mat2str(Matrix2(2,1))), "Sistema utilizado" , "none"); 
    case(4)
      [matrizAproximacion] = aproximacionExponencial(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,5); matrizAproximacion(cantFilas(1,2)+1,4)] ;
      msgbox(disp(matrizAproximacion), "Detalle de calculos:\n\n", "none"); 
      msbox(cstrcat(mat2str(Matrix1(1,1)), " * a + ", mat2str(Matrix1(1,2)), " * B = ", mat2str(Matrix2(1,1)), "\n\n", mat2str(Matrix1(2,1)), " * a + ", mat2str(Matrix1(2,2)), " * B = ", mat2str(Matrix2(2,1))), "Sistema utilizado" , "none");
    case(5)
      [matrizAproximacion] = aproximacionHiperbolica(listX, listY, cantidadDecimales);
      Matrix1 = [cantFilas(1,2),matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1),matrizAproximacion(cantFilas(1,2)+1,3)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,2);matrizAproximacion(cantFilas(1,2)+1,4)];
      msgbox(disp(matrizAproximacion), "Detalle de calculos:\n\n", "none"); 
      msbox(cstrcat(mat2str(Matrix1(1,1)), " * a + ", mat2str(Matrix1(1,2)), " * B = ", mat2str(Matrix2(1,1)), "\n\n", mat2str(Matrix1(2,1)), " * a + ", mat2str(Matrix1(2,2)), " * B = ", mat2str(Matrix2(2,1))), "Sistema utilizado" , "none"); 
  endswitch
endfunction

function graficar(tipoFuncion, listX, listY, cantidadDecimales)
  [funcion]= obtenerFuncion(tipoFuncion,listX,listY, cantidadDecimales);
  ploteoPuntos(listX,listY,funcion);
endfunction

function mostrarFuncionAproximante(tipoFuncion, listX, listY, decimales)
  [funcion]= obtenerFuncion(tipoFuncion,listX,listY, decimales);
  syms x;
  g = symfun(funcion,x);
  msgbox(cstrcat("La funcion es: \n\n",disp(g),"\n\n"));
endfunction

function [funcion] = obtenerFuncion(tipoFuncion,listX,listY, cantidadDecimales)
  cantFilas = size(listX);
  switch(tipoFuncion)
    case(1)
      [matrizAproximacion] = aproximacionLineal(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,4); matrizAproximacion(cantFilas(1,2)+1,2)] ;
      Solucion = Matrix1\Matrix2;
      Solucion(1,1) = redondear(Solucion(1,1),cantidadDecimales);
      Solucion(2,1) = redondear(Solucion(2,1),cantidadDecimales);
      funcion = @(x)x*Solucion(1,1)+Solucion(2,1);
    case(2)
      [matrizAproximacion] = aproximacionParabola(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,5);matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4);cantFilas(1,2), matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,7); matrizAproximacion(cantFilas(1,2)+1,6);matrizAproximacion(cantFilas(1,2)+1,2)];
      Solucion = Matrix1\Matrix2;
      Solucion(1,1) = redondear(Solucion(1,1),cantidadDecimales);
      Solucion(2,1) = redondear(Solucion(2,1),cantidadDecimales);
      funcion = @(x)x.^2*Solucion(1,1)+x*Solucion(2,1)+Solucion(3,1);
    case(3)
      [matrizAproximacion] = aproximacionPotencial(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,3) ; matrizAproximacion(cantFilas(1,2)+1,3), cantFilas(1,2)]
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,6); matrizAproximacion(cantFilas(1,2)+1,5)];
      Solucion = Matrix1\Matrix2;
      sizeMatrix2 = size(Matrix2);
      Solucion(1,1) = redondear(Solucion(1,1),cantidadDecimales);
      Solucion(sizeMatrix2(1,1),1) = redondear(exp(Solucion(sizeMatrix2(1,1),1)),cantidadDecimales);
      Solucion(2,1) = redondear(Solucion(2,1),cantidadDecimales);
      funcion = @(x)Solucion(2,1)*x.^(Solucion(1,1));
    case(4)
      [matrizAproximacion] = aproximacionExponencial(listX, listY, cantidadDecimales);
      Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,5); matrizAproximacion(cantFilas(1,2)+1,4)] ;
      Solucion = Matrix1\Matrix2;
      Solucion(1,1) = redondear(Solucion(1,1),cantidadDecimales);
      Solucion(2,1) = redondear(exp(Solucion(2,1)),cantidadDecimales);
      funcion = @(x)Solucion(2,1)*e.^(Solucion(1,1)*x) ;
    case(5)
      [matrizAproximacion] = aproximacionHiperbolica(listX, listY, cantidadDecimales);
      Matrix1 = [cantFilas(1,2),matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1),matrizAproximacion(cantFilas(1,2)+1,3)];
      Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,2);matrizAproximacion(cantFilas(1,2)+1,4)];
      Solucion = Matrix1\Matrix2;
      Solucion(1,1) = redondear(Solucion(1,1)*(Solucion(2,1).**(-1)),cantidadDecimales);
      Solucion(2,1) = redondear(Solucion(2,1).**(-1),cantidadDecimales);
      funcion = @(x)Solucion(2,1)*((x+Solucion(1,1)).^(-1));
      
  endswitch
endfunction








