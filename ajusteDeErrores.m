function main
  listX = input("Ingrese una lista de valores para x")
  listY = input("Ingrese una lista de valores para y")
  aproximacionParabola(listX, listY)

endfunction

function aproximacionPotencial(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 6)
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i)
    matrizAproximacion(i,2) = listY(1,i)
    matrizAproximacion(i, 3) = reallog(listX(1,i))
    matrizAproximacion(i, 4) = reallog(listX(1,i)).**2
    matrizAproximacion(i,5) = reallog(listY(1,i))
    matrizAproximacion(i,6) = reallog(listX(1,i))*reallog(listY(1,i))
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1))
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2))
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3))
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4))
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5))
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6))  
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,3) ; matrizAproximacion(cantFilas(1,2)+1,3), cantFilas(1,2)]
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,6); matrizAproximacion(cantFilas(1,2)+1,5)]
  Solucion = Matrix1\Matrix2
  sizeMatrix2 = size(Matrix2)
  Solucion(sizeMatrix2(1,1),1) = exp(Solucion(sizeMatrix2(1,1),1))
endfunction

function aproximacionLineal(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4)
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i)
    matrizAproximacion(i,2) = listY(1,i)
    matrizAproximacion(i,3) = listX(1,i).**2
    matrizAproximacion(i,4) = listX(1,i)*listY(1,i)
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1))
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2))
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3))
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4))
  
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1), cantFilas(1,2)]
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,4); matrizAproximacion(cantFilas(1,2)+1,2)] 
  Solucion = Matrix1\Matrix2
endfunction

function aproximacionHiperbolica(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 4)
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i)
    matrizAproximacion(i,2) = listY(1,i)
    matrizAproximacion(i, 3) = listX(1,i).**2
    matrizAproximacion(i, 4) = listX(1,i)*listY(1,i)
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1))
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2))
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3))
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4))
  Matrix1 = [cantFilas(1,2),matrizAproximacion(cantFilas(1,2)+1,1);matrizAproximacion(cantFilas(1,2)+1,1),matrizAproximacion(cantFilas(1,2)+1,3)]
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,2);matrizAproximacion(cantFilas(1,2)+1,4)]
  Solucion = Matrix1\Matrix2
  Solucion(1,1) = Solucion(1,1)*(Solucion(2,1).**(-1))
  Solucion(2,1) = Solucion(2,1).**(-1)
endfunction

function aproximacionParabola(listX, listY)
  cantFilas = size(listX)
  matrizAproximacion = zeros(cantFilas(1,2)+1, 7)
  for i = 1:cantFilas(1,2)
    matrizAproximacion(i,1) = listX(1,i)
    matrizAproximacion(i,2) = listY(1,i)
    matrizAproximacion(i,3) = listX(1,i).**2
    matrizAproximacion(i,4) = listX(1,i).**3
    matrizAproximacion(i,5) = listX(1,i).**4
    matrizAproximacion(i,6) = listX(1,i)*listY(1,i)
    matrizAproximacion(i,7) = listY(1,i)*listX(1,i).**2
  endfor
  matrizAproximacion(cantFilas(1,2)+1,1) = sum(matrizAproximacion(:,1))
  matrizAproximacion(cantFilas(1,2)+1,2) = sum(matrizAproximacion(:,2))
  matrizAproximacion(cantFilas(1,2)+1,3) = sum(matrizAproximacion(:,3))
  matrizAproximacion(cantFilas(1,2)+1,4) = sum(matrizAproximacion(:,4))
  matrizAproximacion(cantFilas(1,2)+1,5) = sum(matrizAproximacion(:,5))
  matrizAproximacion(cantFilas(1,2)+1,6) = sum(matrizAproximacion(:,6))
  matrizAproximacion(cantFilas(1,2)+1,7) = sum(matrizAproximacion(:,7))
  Matrix1 = [matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4), matrizAproximacion(cantFilas(1,2)+1,5);matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3), matrizAproximacion(cantFilas(1,2)+1,4);cantFilas(1,2), matrizAproximacion(cantFilas(1,2)+1,1), matrizAproximacion(cantFilas(1,2)+1,3)]
  Matrix2 = [matrizAproximacion(cantFilas(1,2)+1,7); matrizAproximacion(cantFilas(1,2)+1,6);matrizAproximacion(cantFilas(1,2)+1,2)]
  Solucion = Matrix1\Matrix2
endfunction