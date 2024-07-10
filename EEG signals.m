clc
clear
close all

%% PARTE PRIMERA %% 
% -	Realice la carga de la base de datos propuesta para esta actividad utilizando la columna 3, en Matlab.
% -	Realice segmentaciones de datos cada 6000 muestras hasta el fin del
% registro con una frecuencia de muestreo de 2048 Hz
% -	Compare los Espectros de Densidad de Potencia de 5 intervalos equiespaciados de elementos muestreados con una representación gráfica logarítmica.

%Primero se cargan los datos del archivo Trial1
MyCell = readcell('Trial1.csv'); %lee el archivo Trial1 de los datos del EEG y lo almacena en MyCell
x=cell2mat(MyCell(:,3)); %Almacena la tercera columna de MyCell en 'x'

%Se define la longitud del segmento y la frecuencia de muestreo
L=6000; %Longitud del segmento: 6000 muestras
fs=2048; %Frecuencia de muestreo: 2048Hz

%Se calcula el numero de segmentos
Segmentos=round(length(x)/L); %Agrupa los datos en intervalos de 6000, redondeando la longitud de la tercera columna x, dividida por 6000 muestras

%Se inicializa la matriz B para almacenar segmentos
B=zeros(L,Segmentos); %Se crea una matriz B de ceros de las dimensiones de la tercera columna de MyCell y el numero de segmentos calculado antes

%Segmentar los datos
for i=1:Segmentos
    B(:,i)=x(L*(i-1)+1:L*i);
end
%Se utiliza un bucle for que divide los datos de x en segmentos de 6000 muestras y los almacena en la matriz B inicializada antes

% Selección de 5 segmentos para análisis de intervalos (equidistantes cada 47)
seg1  =   B(:,1);   % Numero de segmento 1 
seg2  =   B(:,47);  %  Numero de segmento 47
seg3  =   B(:,94);  % Numero de segmento 94 (47+47)
seg4  =   B(:,141); % Numero de segmento 141  (94+47)
seg5  =   B(:,188); % Numero de segmento 188 (141+47)
segseleccionados=[seg1, seg2, seg3, seg4, seg5];
sc=size(segseleccionados,2);

f=1; %Inicializa la variable f, variable de las figuras

%Se va a realizar los espectros de potencia en eje logaritmico de cada segmento en DC

%Para comparar los espectros de potencia necesitamos hacer la transformada de Fourier de los 5 segmentos, para ello utilizamos un bucle for
for i=1:sc  
    
    Z=segseleccionados(:,i); %Coje el i-esimo segmento de la matriz segseleccionados y lo asigna a la variable Z
    ZT=transpose(Z); %Hace la transpuesta del vector Z
        
    Y=fft(ZT); %Transformada rapida de Fourier de la matriz transpuesta
    Pt = abs(Y); %Se obtiene el valor absoluto, el modulo de la transformada
    
    %Se crean los ejes de frecuencia y potencia
    f_axis = 0:(fs/L):(fs/2-fs/L); %genera valores equiespaciados entre 0 y la mitad de la frecuencia de muestreo con una separacion de fs/L
    P = Pt(1:L/2); % %valor absoluto, para calcular la magnitud del espectro de potencia

    %Se reperentan los espectros 
    figure(f) %Crea una nueva figura
    semilogx(f_axis,P,'b-')  %Realiza la gráfica en escala semilogarítmica, color azul y linea continua
    xlabel('Frecuencia (Hz)') %Titulo del eje x
    ylabel('Magnitud |P(f)|') %Titulo del eje y
    title('Espectro de Potencia (logarítmico) para el segmento seleccionado', i ) %Titulo de la grafica
    grid on %Cuadricula

    f=f+1; %Aumenta el numero de la figura

end


%% PARTE SEGUNDA %% 
% -	Elimine el nivel DC de la señal (a cada uno de los segmentos anteriores elimine el promedio de la señal).
% -	Halle la Desviación Estándar de la Señal en los segmentos señalados. 
% -	Compare dichas desviaciones. 
% -	Realice los Espectros de Densidad de Potencia y compárelos con los de la primera actividad, verifique ¿Qué cambios existen? Si no observa cambios significativos pruebe a realizar dicha comparativa en representación lineal en vez de logarítmica.

media=zeros(size(segseleccionados,2)); %Inicializa un vector media de tamaño igual al numero de segmentos seleccinados con ceros para almacenar las medias de cada segmento
desvstd=zeros(size(segseleccionados,2)); %Inicializa un vector desvstd de tamaño igual al numero de segmentos seleccionados con ceros para almacenar las desviaciones estandar de cada segmento
     
%Primero, se realizan los Espetros de Potencia en ejes logaritmicos de cada segmento sin DC 
for i=1:sc  

    Z=segseleccionados(:,i); %Extrae el i-esimo segmento de la matriz segseleccionados y lo asigna a la variable z

    %Eliminamos el nivel DC, para ello buscamos el promedio de los segmentos anteriores y lo eliminamos
    media(i)=mean(Z); %calcula la media de los valores de Z y la almacena en media. Esto representa el nivel DC de la señal
    ZAC=Z-media(i); %Resta la media calculada antes al vector Z, par eliminar el nivel DC de la señal
    
    desvstd(i)=std(ZAC); %calcula la desviacion estandar de la señal
    
    %Realizamos el Espectro de potencia sin DC de la señal 
    ZTAC=transpose(ZAC); %Realizamos la transpuesta del segmento sin DC   
    YAC=fft(ZTAC); %Hacemos la transformada Rapida de Fourier del segmento sin DC transpueto
    PtAC = abs(YAC); %Calcula el valor absoluto

    %Se crean los ejes de frecuencia y potencia
    f_axis = 0:(fs/L):(fs/2-fs/L); %genera valores equiespaciados entre 0 y la mitad de la frecuencia de muestreo con una separacion de fs/L
    PAC = PtAC(1:L/2); %valor absoluto, para calcular la magnitud del espectro de potencia

    %Se representan los espectros sin el nivel DC en escala semilogaritmica
    figure(f) %Crea una nueva figura
    semilogx(f_axis,PAC,'r-') %Realiza la gráfica en escala semilogarítmica de la magnitud del espectro de potencia sin DC en funcion de la frecuencia, color rojo y linea continua
    xlabel('Frecuencia (Hz)') %Titulo del eje x
    ylabel('Magnitud |P(f)|') %Titulo del eje y
    title('Espectro de Potencia sin DC (logarítmico) para el segmento seleccionado', i ) %Titulo de la grafica
    grid on %Cuadricula

    f=f+1;

end
%Como no se observa diferencia con las escalas logaritmicas, se hacen en lineal

%Para realizar el espectro de potencia lineal para cada segmento sin
%haberles quitado la DC, se repite el procedimiento anterior, calculando
%la transpusta, tranformada de fourier, el valor absoluto, los ejes, y se
%representa de forma lineal

for i=1:sc  

    Z=segseleccionados(:,i);
    ZT=transpose(Z);

    Y=fft(ZT);

    Pt = abs(Y);

    f_axis = 0:(fs/L):(fs/2-fs/L);
    P = Pt(1:L/2);

    figure(f)
    plot(f_axis,P,'b-') %Realiza un grafico lineal de la magnitud del Espectro de potencia en funcion de la frecuencia, la linea es azul y continua
    xlabel('Frecuencia (Hz)')
    ylabel('Magnitud |P(f)|')
    title('Espectro de Potencia (lineal) para el segmento seleccionado', i )
    grid on

    f=f+1;

end

%Por ultimo, se realiza el espectro de potencia para cada segmento
%sin DC en una representacion lineal
for i=1:sc  

    Z=segseleccionados(:,i);

    media(i)=mean(Z);
    ZAC=Z-media(i);
    desvstd(i)=std(ZAC);

    ZTAC=transpose(ZAC);
        
    YAC=fft(ZTAC);
    PtAC = abs(YAC);

    f_axis = 0:(fs/L):(fs/2-fs/L);
    PAC = PtAC(1:L/2);

    figure(f)
    plot(f_axis,PAC,'r-')

    xlabel('Frecuencia (Hz)')
    ylabel('Magnitud |P(f)|')
    title('Espectro de Potencia sin DC (lineal) para el segmento seleccionado', i )
    grid on

    f=f+1;

end


%% PARTE TERCERA %% 
% -	Implemente filtrado digital con la finalidad de hallar los ritmos Delta, Theta, Alfa, Beta y Gama. 
% -	Realice los Espectros de Densidad de Potencia de cada ritmo, comente
% los resultados de ello, ¿coincide con lo esperado?

%Se establecen variables igual a cero para analizar los valors minimos y
%maximos de las magnitudes del tiempo y frecuencia
min_magnitud_t = 0;
max_magnitud_t = 0;
max_magnitud_f = 0;

%Se utiliza un bucle for que recorre cada segmento Z y le resta la media
%para eliminar la corriente continua
for i = 1:sc
    Z = segseleccionados(:, i);
    Z = Z - mean(Z);

%Se definen los filtros de banda. Se especifican las frecuencias de corte
%para cada ritmo cerebral
    % Filtros característicos para EEG
    delta_banda = [0.5 4];    % Frecuencias cerebral Delta
    theta_banda = [4 8];    %  Frecuencias cerebral Theta
    alpha_banda = [8 13];    %  Frecuencias cerebral Alpha
    beta_banda  = [13 30];    %  Frecuencias cerebral Beta
    gamma_banda = [30 100];    %  Frecuencias cerebral Gamma

%Se filtran las señales. Se utiliza un filtro digital de tipo bandpass para 
% obtener la señal en las frecuencias de cada ritmo
    delta_filtrada = bandpass(Z, delta_banda, fs);
    theta_filtrada = bandpass(Z, theta_banda, fs);
    alpha_filtrada = bandpass(Z, alpha_banda, fs);
    beta_filtrada  = bandpass(Z, beta_banda, fs);
    gamma_filtrada = bandpass(Z, gamma_banda, fs);
    
%Se calculan las transformadas de Fourier de las señales filtradas en cada
%banda y se utiliza el valor absoluto para obtener la magnitud
    delta_Pt = abs(fft(delta_filtrada));
    theta_Pt = abs(fft(theta_filtrada));
    alpha_Pt = abs(fft(alpha_filtrada));
    beta_Pt  = abs(fft(beta_filtrada));
    gamma_Pt = abs(fft(gamma_filtrada));

%Se calculan los valores  minimo y maximo de amplitud entre todas las 
% señales filtradas en el dominio del tiempo.
% Y el maximo para el dominiuo de la frecuencia
    min_magnitud_banda_t = min([min(delta_filtrada), min(theta_filtrada), min(alpha_filtrada), min(beta_filtrada), min(gamma_filtrada)]);
    max_magnitud_banda_t = max([max(delta_filtrada), max(theta_filtrada), max(alpha_filtrada), max(beta_filtrada), max(gamma_filtrada)]);
    max_magnitud_banda_f = max([max(delta_Pt), max(theta_Pt), max(alpha_Pt), max(beta_Pt), max(gamma_Pt)]);

%Para que los espectros esten escalados de manera adecuada: 
%Se coje el valor minimo entre el valor min_magnitud_t y el calculado antes
%min_magnitud_banda_t, se hace lo mismo para los maximos en el dominio del
%tiempo y de la frecuencia
    min_magnitud_t = min(min_magnitud_t, min_magnitud_banda_t);
    max_magnitud_t = max(max_magnitud_t, max_magnitud_banda_t);
    max_magnitud_f = max(max_magnitud_f, max_magnitud_banda_f);

    %Se crean los ejes de frecuencia
    f_axis = 0:(fs / L):(fs / 2 - fs / L);

    %Se crea una nueva figura 
    figure(f)

    %A continuacion se van a crear figuras con subgraficas.
    
    %Crea la primera subgráfica: la señal original
    subplot(7, 5, 1:5); %Para organizar las graficas, se divide la figura 
        % en 7 filas y 5 columnas, y esta gráfica ocupará todas las columnas
    plot(1:length(Z), Z, 'b-') %Hace la gráfica de la señal original, con 
        % rango 1 a la longitud de la señal. Linea azul continua
    title(['Señal Original - Segmento ', num2str(i)]) %Añade el titulo de 
        % la gráfica y con num2str el numero de segmento
    ylabel('Amplitud') %Titulo del eje y: Amplitud
    xlabel('Muestras') %Titulo del eje x: Muestras
    legend('Señal original') %Añade una leyenda
    grid on %Cuadricula

    %Crea la segunda subgrafica: la señal filtrada en la banda de frecuencia
    % Delta
    subplot(7, 5, 6:10); %
    plot(1:length(delta_filtrada), delta_filtrada, 'g-')
    title(['Delta - Segmento ', num2str(i)])
    ylabel('Amplitud')
    ylim([min_magnitud_t max_magnitud_t]) %utiliza el limites de eje 
        % calculados antes, para facilitar la comparación
    xlabel('Muestras')
    legend('Delta')
    grid on

    %Crea la tercera subgráfica: la señal filtrada en la banda de frecuencia 
    % Theta
    subplot(7, 5, 11:15);
    plot(1:length(theta_filtrada), theta_filtrada, 'r-')
    title(['Theta - Segmento ', num2str(i)])
    ylabel('Amplitud')
    ylim([min_magnitud_t max_magnitud_t])
    xlabel('Muestras')
    legend('Theta')
    grid on

    %Crea la cuarta subgráfica: la señal filtrada en la banda de frecuencia 
    % Alpha
    subplot(7, 5, 16:20);
    plot(1:length(alpha_filtrada), alpha_filtrada, 'c-')
    title(['Alpha - Segmento ', num2str(i)])
    ylabel('Amplitud')
    ylim([min_magnitud_t max_magnitud_t])
    xlabel('Muestras')
    legend('Alpha')
    grid on

    %Crea la quinta subgráfica: la señal filtrada en la banda de frecuencia
    % Beta
    subplot(7, 5, 21:25);
    plot(1:length(beta_filtrada), beta_filtrada, 'm-')
    title(['Beta - Segmento ', num2str(i)])
    ylabel('Amplitud')
    ylim([min_magnitud_t max_magnitud_t])
    xlabel('Muestras')
    legend('Beta')
    grid on

    %Crea la sexta subgráfica: la señal filtrada en la banda de frecuencia 
    % Gamma
    subplot(7, 5, 26:30);
    plot( 1:length(gamma_filtrada), gamma_filtrada, 'y-')
    title(['Gamma - Segmento ', num2str(i)])
    ylabel('Amplitud')
    ylim([min_magnitud_t max_magnitud_t])
    xlabel('Muestras')
    legend('Gamma')
    grid on

    %Ahora se van a crear las subgraficas de los espectros de densidad 
    % de Potencia de cada banda de frecuencia filtrada
    %Espectro de densidad de Potencia de la banda de frecuencia Delta
    subplot(7, 5, 31);
    semilogx(f_axis, delta_Pt(1:L/2), 'b-') %Utiliza la escala 
    % semilogaritmica,
    % el eje x son las frecuencias, y el eje y las magnitudes
    title(['Delta - Segmento ', num2str(i)])
    ylabel('Magnitud |P(f)|')
    ylim([0 max_magnitud_f]) %Límites del eje y, las frecuencias, 
    % el valor minimo será 0 y el máximo el calculado antes
    xlabel('Frecuencia (Hz)')
    grid on
 
    %Espectro de densidad de Potencia de la banda de frecuencia Theta
    subplot(7, 5, 32);
    semilogx(f_axis, theta_Pt(1:L/2), 'b-')
    title(['Theta - Segmento ', num2str(i)])
    ylabel('Magnitud |P(f)|')
    ylim([0 max_magnitud_f])
    xlabel('Frecuencia (Hz)')
    grid on
    
    %Espectro de densidad de Potencia de la banda de frecuencia Alpha
    subplot(7, 5, 33);
    semilogx(f_axis, alpha_Pt(1:L/2), 'b-')
    title(['Alpha - Segmento ', num2str(i)])
    ylabel('Magnitud |P(f)|')
    ylim([0 max_magnitud_f])
    xlabel('Frecuencia (Hz)')
    grid on

    %Espectro de densidad de Potencia de la banda de frecuencia Beta
    subplot(7, 5, 34);
    semilogx(f_axis, beta_Pt(1:L/2), 'b-')
    title(['Beta - Segmento ', num2str(i)])
    ylabel('Magnitud |P(f)|')
    ylim([0 max_magnitud_f])
    xlabel('Frecuencia (Hz)')
    grid on

    %Espectro de densidad de Potencia de la banda de frecuencia Gamma
    subplot(7, 5, 35);
    semilogx(f_axis, gamma_Pt(1:L/2), 'b-')
    title(['Gamma - Segmento ', num2str(i)])
    ylabel('Magnitud |P(f)|')
    ylim([0 max_magnitud_f])
    xlabel('Frecuencia (Hz)')
    grid on

    f = f + 1;
end



