%% Rensa och ladda in data.
clear all;
close all;


%% L�s in bilderna
% OBS! Bilderna m�ste vara i ordning, s� bilderna m�ste heta Img01, Img02
% osv.. ist�llet f�r Img1, Img2.. osv.
directory = dir('N*.jpg');
num = numel(directory);
images = cell(1, num);
imsmall = cell(1, num);
imbw = cell(1, num);

for n = 1:num
    images{n} = im2double(imread(directory(n).name));
    imsmall{n} = imresize(images{n}, 0.5);
    imbw{n} = rgb2gray(imsmall{n});
end

% Storlek p� bilderna, X-bredd, Y-h�jd.
[Y, X] = size(imbw{1});


%% Visa bilderna
montage([imsmall{1}, imsmall{2}, imsmall{3}]);


%% Hitta kanterna
imedges = cell(1, num);

% Antal punkter att hitta.
npoints = 200;

for n = 1:num
    % Sortera alla h�rn i bilderna f�r att kunna v�lja de punkter i b�da
    % kanterna av bilderna. Dessa kanter hittas med Harris-detekorn
    % OBS! "corners" fungerar troligtvis bara i R2013.
    sorted_corners = sortrows(corner(imbw{n}, 'Harris'), 1);
    imedges{n} = [sorted_corners(1:npoints, 1) sorted_corners(1:npoints, 2) ;
               sorted_corners(end-npoints+1:end, 1) sorted_corners(end-npoints+1:end, 2)];
end


%% Ett f�rs�k att hitta intressanta punkter utan en automagisk Harrisd.
% -- FUNGERAR INTE SOM DET VAR T�NKT, K�R EJ DENNA. --
pfilt = [1 1 1;
         1 -8 1;
         1 1 1];

for n = 1:num
    imbw{n} = rgb2gray(imsmall{n});
    randomint = conv2(imbw{n}, pfilt);
    imedges{n} = randomint;
end

imshow(imedges{1}, []);

[x, y] = find(imedges{1} == max(max(imedges{1})))


%% Visa dessa intressanta punkter i samtliga bilder.
close all

for n = 1:num
    figure, imshow(imbw{n}, []);
    hold on, plot(imedges{n}(:,1), imedges{n}(:,2), '.');
end


%% H�r under �r ett annat f�rs�k till att hitta intressanta punkter
%  Denna metod anv�nder korrelation.
%  F�r n�rvarande anv�nder denna metoden punkter fr�n Harris detektorn fr�n
%  tidigare f�r att hitta intressanta punkter att anv�nda korrelation
%  mellan bilderna.
%
%
%
%
%
%
%
%
%
%
%
clc;

% Endast 15% av bilderna �r intressanta att titta p� (�verlappning).
impercent = 0.85;

% Hur stora korrelationsbilderna ska vara, ex. 15 blir (15*2+1)^2 pixlar
corrimspace = 15;

% Allokerna minne f�r de gamla koordinaterna.
oldcoords = zeros(5, 2);

% Skapa en pool av intressanta punkter.
poi = imedges{1}(npoints+1:end, :);

% Initiera iterator f�r while-loopen nedan.
iter = 0;

while iter < 5
    iter = iter + 1;

    % Skaffa fram en godtycklig koordinat.
    randomint = randi(size(poi, 1));
    randomcoord = poi(randomint, :);

    % Kolla om vi kan skapa en korrelationsbild fr�n denna punkten.
    if (randomcoord(1) > X - corrimspace || randomcoord(2) > Y - corrimspace || randomcoord(1) < floor(X*impercent))
        iter = iter - 1;
        continue;
    end

    % Spara denna godtyckliga koordinat i en ny vektor.
    oldcoords(iter, 1) = poi(randomint, 1);
    oldcoords(iter, 2) = poi(randomint, 2);

    % Ta bort denna koordinat fr�n poolen.
    poi(randomint, :) = [];
end


%% Visa f�rsta bilden med intressanta punkter.
close all

figure, imshow(imbw{1}, []);
hold on, plot(oldcoords(:,1), oldcoords(:,2), '.');


%% Skapa korrelationsbilder f�r de intressanta punkterna.
crop = cell(1, size(oldcoords, 1));
for n = 1:size(oldcoords, 1)
    crop{n} = imbw{1}(oldcoords(n, 2)-corrimspace:oldcoords(n, 2)+corrimspace, oldcoords(n, 1)-corrimspace:oldcoords(n, 1)+corrimspace);
    figure, imshow(crop{n});
end


%% Hitta punkterna i bild 2 med korrelation.
newcoords = zeros(5, 2);

% Kolla bara p� en delm�ngd av bilden.
small_image2 = imbw{2};
small_image2(:, floor(X*(1-impercent)):end) = [];

for n = 1:5
    temp = normxcorr2(crop{n}, small_image2);

    [y, x] = find(temp == max(temp(:)));
    newcoords(n, 1) = x;
    newcoords(n, 2) = y;
    figure, imshow(temp, []);
    hold on, plot(x, y, 'r.');
end


%% Sannlingens �gonblick, visa intressanta punkter i bild tv�, st�mmer dem?
figure, imshow(imbw{2}, [])
hold on, plot(newcoords(:, 1), newcoords(:, 2), 'r.');


%% Nu ska vi testa om koordinaterna var n�got att ha genom att g�ra en
%  m�ngd iterationer
%
%
%
%
%
%
%
%
%
%
%
clc;

% Som tidigare beh�vs endast 15%.
impercent = 0.85;

% Vi ska ta fram s�h�r m�nga transformationsmatriser.
passes = 1000;

% H�r sparas alla transformationer.
xes = cell(1, passes);

% Skapa tv� pooler, en f�r varje bild.
poi = imedges{1}(npoints+1:end, :);
poi2 = imedges{2}(1:npoints, :);

for k = 1:passes

    % Spara koordinaterna i denna bilden, ny och fr�sch f�r varje iter.
    oldcoords = zeros(5, 2);
    oldcoords2 = zeros(5, 2);

    % Spara poolen med intressanta punkter, ny och fr�sch f�r varje iter.
    poi_temp = poi;
    poi_temp2 = poi2;

    iter = 0;

    while iter < 5

        iter = iter + 1;

        % Hitta slumpm�ssigt ett koordinatpar.
        randomint = randi(size(poi_temp, 1));
        randomcoord = poi_temp(randomint, :);

        randomint2 = randi(size(poi_temp2, 1));
        randomcoord2 = poi_temp2(randomint2, :);

        % Spara dessa.
        oldcoords(iter, 1) = poi_temp(randomint, 1);
        oldcoords(iter, 2) = poi_temp(randomint, 2);

        oldcoords2(iter, 1) = poi_temp2(randomint2, 1);
        oldcoords2(iter, 2) = poi_temp2(randomint2, 2);

        % Ta bort dem fr�n poolen.
        poi_temp(randomint, :) = [];
        poi_temp2(randomint, :) = [];

    end

    % Skapa sedan en transformation f�r varje koordinatpar.
	xes{k} = [oldcoords ones(5, 1)] \ [oldcoords2  ones(5, 1)];

end


%% Nu har vi en lista med transformationer som st�mmer med bilderna p�
%  olika bra s�tt, nu ska vi kolla kvalit�n p� dessa med v�ra
%  kontrollpunkter, dvs. alla andra punkter som vi tagit ut som h�rn
%  tidigare.

% Vi ska kolla s�h�r m�nga g�nger.
passes = 1000;

% H�r sparas det om det �r en bra transform.
wasitgoodindex = zeros(1, passes);

for k = 1:passes

    oldcoords = zeros(10, 2);
    oldcoords2 = zeros(10, 2);

    poi_temp = poi;
    poi_temp2 = poi2;

    for iter = 1:20

        % Get a random coordinate
        randomint = randi(size(poi_temp, 1));
        randomcoord = poi_temp(randomint, :);

        randomint2 = randi(size(poi_temp2, 1));
        randomcoord2 = poi_temp2(randomint2, :);

        % Not the points is X, Y as usual. Jeeez.
        oldcoords(iter, 1) = poi_temp(randomint, 1);
        oldcoords(iter, 2) = poi_temp(randomint, 2);

        oldcoords2(iter, 1) = poi_temp2(randomint2, 1);
        oldcoords2(iter, 2) = poi_temp2(randomint2, 2);

        % Remove the picked poi from the poi pool.
        poi_temp(randomint, :) = [];
        poi_temp2(randomint, :) = [];

    end
    

    % Transform our P1 coords to ~P2 coords
    tempcoords = [oldcoords ones(20, 1)] * xes{k};
    
    eukdist = 0;
    for n = 1:20
        eukdist = eukdist + sqrt((tempcoords(n, 1)-oldcoords2(n, 1))^2 + (tempcoords(n, 2)-oldcoords2(n, 2))^2);
    end
    
    wasitgoodindex(k) = eukdist;

end



rightindex = find(wasitgoodindex == min(wasitgoodindex(:)));

theX = xes{rightindex};

newPic = 0;

for i = 1:X
    for j = 1:Y
        uvw1 = [i, j, 1] * theX;
        u = round(uvw1(1) / uvw1(3));
        v = round(uvw1(2) / uvw1(3));
        if (u > 0 && u <= X && v > 0 && v <= Y)
            newPic(u, v) = imbw{1}(j, i);
            filled(u, v) = 1;
        end
    end
end

figure, imshow(newPic, []);




%% H�r gav vi upp lite och k�r manuella punkter med ginput om s� �nskas.
%  Det finna f�rdiga koordinater l�ngre ner i koden, s� denna beh�vs ej
%  k�ras.
%
%
%
%
%
%
%

% Antal punkter som ska plockas.
numcoords = 6;

% Plocka punkter med ginput.
% Fr�n Bild 1 & Bild 2.
im1_poi = zeros(1, 2);
im2_poi = zeros(1, 2);

for i = 1:numcoords
    
    imshow(imbw{1});
    [x_temp, y_temp] = ginput;
    
    im1_poi = [im1_poi; [x_temp y_temp]];

    imshow(imbw{2});
    
    [x_temp, y_temp] = ginput;
    
    im2_poi = [im2_poi; [x_temp y_temp]];

end

close all;

im1_poi(1,:) = [];
im2_poi(1,:) = [];

% Krydda koordinaterna med homogena kompisar.
im1_poi = [im1_poi ones(numcoords,1)];
im2_poi = [im2_poi ones(numcoords,1)];


%% Pss fast med Bild 2 -> Bild 3
im23_poi = zeros(1, 2);
im3_poi = zeros(1, 2);

for i = 1:numcoords
    
    imshow(imbw{2});
    [x_temp, y_temp] = ginput;
    
    im23_poi = [im23_poi; [x_temp y_temp]];

    imshow(imbw{3});
    
    [x_temp, y_temp] = ginput;
    
    im3_poi = [im3_poi; [x_temp y_temp]];

end

close all;

im23_poi(1,:) = [];
im3_poi(1,:) = [];

im23_poi = [im23_poi ones(size(im23_poi, 1),1)];
im3_poi = [im3_poi ones(size(im3_poi, 1),1)];


%% K�r lite goa koordinater annars.
im1_poi = [1227,25,1;1110,166,1;1084,356,1;1204,724,1;1191,110,1;1250,488,1];
im2_poi = [153,64,1;39,181,1;11,375,1;129,748,1;120,139,1;171,513,1];
im23_poi = [1124,258,1;1241,281,1;1179,365,1;1143,518,1;1215,872,1;1276,402,1];
im3_poi = [44,273,1;165,310,1;101,388,1;62,542,1;119,891,1;193,429,1];


%% Om man vill kontrollera transformationen s�... Bild 1 -> Bild 2
trans = im1_poi\im2_poi;

plot(im1_poi(:, 1), im1_poi(:, 2), 'ro');

hold on
plot(im2_poi(:, 1), im2_poi(:, 2), 'go');

Y_maybe = im1_poi*trans;
plot(Y_maybe(:, 1), Y_maybe(:, 2), 'bo');


%% Om man vill kontrollera transformationen s�... Bild 3 -> Bild 2
trans2 = im3_poi\im23_poi;

plot(im23_poi(:, 1), im23_poi(:, 2), 'ro');

hold on
plot(im3_poi(:, 1), im3_poi(:, 2), 'go');

Y_maybe = im3_poi*trans2;
plot(Y_maybe(:, 1), Y_maybe(:, 2), 'bo');



%% H�r skall bilderna transformeras och s�ttas ihop till en!

% Padda mittenbilden s� att man kan f� in de andra till resultatbilden.
res = padarray(imbw{2}, [0 X]);

[resy, resx] = size(res);

% Spara en matris s� att man kan h�lla koll p� vilka pixlar som blivit
% fyllda. Detta f�r att interpolera lite f�rg �ver dessa.
filled = padarray(ones(Y, X), [0 X]);

% G� genom hela resultatbilden.
for x = 1:X
    for y = 1:Y

        % Transformera koordinaterna f�r den f�rsta bilden
        uvw1 = [x-X, y, 1] * trans;
        u = round(uvw1(1) / uvw1(3))-40;
        v = round(uvw1(2) / uvw1(3))+115;

        % H�r fylls bilden efter villkor.
        if (u > -X*2 && u <= X && v > 0 && v <= Y)
            res(v, u+X*2) = imbw{1}(y, x);
            filled(v, u+X*2) = 1;
        end
        

        % Pss fast bild 3 -> bild 2
        uvw1 = [x+X, y, 1] * trans2;
        u = round(uvw1(1) / uvw1(3));
        v = round(uvw1(2) / uvw1(3))+90;

        if (u > 2*X && u <= 3*X && v > 0 && v <= Y)
            res(v, u) = imbw{3}(y, x);
            filled(v, u) = 1;
        end
    end
end

imshow(res);


%% N�gon slags interpolering f�r att slippa staketm�sntret.
for x = 1:resx
    for y = 1:resy
        if (filled(y,x) == 0 && x > 1 && x < resx && y > 1 && y < resy)
            nearest = res(y-1:y+1, x-1:x+1);
            nearest_filled = filled(y-1:y+1, x-1:x+1);

            res(y, x) = sum(sum(nearest .* nearest_filled)) / sum(nearest_filled(:));
        end
    end
end

res = res / max(res(:));

figure, imshow(res);




