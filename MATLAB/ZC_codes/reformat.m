function [long,datos]=reformat();
% Programa para leer datos de archivos .dx, .jac ó .txt
%
% La sintaxis es:
%
%	 [long,datos]=reformat;
%
% donde:
%	long  es la variable que contiene las longitudes de onda
%	datos es la variable que contiene los datos espectroscópicos
%
% Versión 1.2 (lee datos de absorción de los archivos *.jac)


disp('Los archivos han de tener un nombre del tipo xxxx01.dx; xxxxx01.jac ó xxxx01.txt')
filein=input('Cadena común de la serie de ficheros de datos (xxxx) ','s');

exten=input('Cadena de la extensión (dx/jac/txt) ','s');
if exten(1:2)=='ja'
	disp('***********************************************************************')
	disp('Este programa está pensado para los nuevos archivos .jac, en los cuales')
	disp('  los datos empiezan en la línea 30.')
	disp('  Si no es éste el caso, se ha de modificar el programa')
	disp('***********************************************************************')
end

nesp=input('Número de ficheros en la serie ');

for j=1:nesp
	if j<10
		filein1=[filein,'0',num2str(j),'.',exten];
	elseif j>=10
		filein1=[filein,num2str(j),'.',exten];
	end
	fid1=fopen(filein1,'r');
	disp(['Estoy abriendo el fichero ',num2str(filein1)])
	n=0;
	m=0;
	while n~=2
		string=fgetl(fid1);
		char=string(1);
		if char =='#' | char=='$'
			m=m+1;
		elseif char~='#'
			n=2;
		end
	end
	frewind(fid1);

	% Excepción en los nuevos ficheros 'jac'. Modificar si procede.
	if exten(1:2)=='ja'
		m=29;
	end

	for i =1:m
		fgetl(fid1);
	end
	data=fscanf(fid1,'%g %g',[2,inf]);
	if (exten(1:2)=='dx' | exten(1:2)=='tx')
		long=data(1,:)';data=data(2,:);
	elseif (exten(1:2)=='ja')
		data=[data(1,:),data(2,:)];
	end
	datos=[datos;data];
	fclose(fid1);
	disp(['Estoy cerrando el fichero ',num2str(filein1)])
	disp('---')
	clear filein1
end

if exten(1:2)=='ja'
	datos=fliplr(datos);
end


% Subrutina extraida del programa matspf.m (escrito por el Gran Tauler)
formato=input('¿Quieres guardar la matriz de datos en formato SQUAD? (s/n) ','s');
if (formato=='s' | formato=='S')
	filename=input('Nombre del fichero SQUAD ','s');
	[n,m]=size(datos);
	fid=fopen(filename,'w');
	for i=1:n;
		for j=1:8:m;
			jfin=j+7;
			if jfin>m;
				jfin=m;
			end
			count1=fprintf(fid,'%10.4f',datos(i,j:jfin)');
			count2=fprintf(fid,'\r\n');
		end
		count3=fprintf(fid,'\r\n');
	end
	status=fclose(fid);
elseif formato ~='s'|'S'
%	break
end

% Esto sí que es una chorrada (escrita por mí, obviamente)
x=round(4*rand)+1;
if x==0
	disp('Después de todo... mañana será otro dia')
elseif x==1
	disp('Ojos que no ven... batacazo que te arreas')
elseif x==2
	disp('Tiran más dos ... que dos carretas')
elseif x==3
	disp('Al fin y al cabo... cada uno baja las escaleras como quiere')
elseif x==4
	disp('Copyright: Grupo de Quimiometria, Diciembre de 1999')
end






