function [number, exponent] = scientific_notation(data)
% Parses data into scientific notation format useful for texts

exponent = floor(log10(abs(data)));
number = sign(data)*round(abs(data/10^exponent));
if number == 10
    number = number/10;
    exponent = exponent + 1;
end