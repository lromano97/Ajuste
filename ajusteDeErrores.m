function main
  listX = input("Ingrese una lista de valores para x")
  listY = input("Ingrese una lista de valores para y")
  aproximacionPotencial(listX, listY)

endfunction

function aproximacionPotencial(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2), 4)
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i, 1) = reallog(listX(1,i))
    matrizAproximacion(i, 2) = reallog(listX(1,i)).**2
    matrizAproximacion(i,3) = reallog(listY(1,i))
    matrizAproximacion(i,4) = reallog(listX(1,i))*reallog(listY(i,1))
  endfor
  sumatoriaLNX = sum(matrizAproximacion(:,1))
  sumatoriaLNXCuadrado = sum(matrizAproximacion(:,2));
  sumatoriaLNY = sum(matrizAproximacion(:,3));
  sumatoriaLNXLNY = sum(matrizAproximacion(:,4));
endfunction


function aproximacionHiperbolica(listX, listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,1), 4);


endfunction

function aproximacionLineal(listX, listY)
  cantFilas = size(listX);
  matrizAproximacion = zeros(cantFilas(1,1), 2);
  for i = 1:cantFilas
    matrizAproximacion(i,1) = listX(i,1).**2;
    matrizAproximacion(i,2) = listX(i,1)*listY(i,1);
  endfor
  sumatoriaXCuadrado = sum(matrizAproximacion(:,1));
  sumatoriaXY = sum(matrizAproximacion(:,2));
  sumatoriaX = sum(listX(:,1));
  sumatoriaY = sum(listY(:,1));
endfunction

function aproximacionParabola(listX, listY)
  cantFilas = size(listX);


endfunction
