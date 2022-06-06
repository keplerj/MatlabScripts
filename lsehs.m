function [mag, ang] = lsehs(f, t, X)
% LSEHS, or least squares extraction of harmonic signals, performs least 
%       squares approximation to given signal(s) X. f is the frequencies to 
%       be approximated for and t the time vector. If X is one signal the 
%       approximation is done in multi-tone fashion. If X is a matrix each 
%       column is assumed to correspond to a different frequency.
%
%       (c) Georgios Printezis 2022

    assert(height(f) == 1)

    assert(width(t) == 1)
    
    if width(X) == 1
        Phi = [cos(2 * pi * f .* t), sin(2 * pi * f .* t), ones(height(t), 1)];
        
        res = Phi \ X;
        
        a = res(1 : width(f));

        b = res(width(f) + 1 : end - 1);

        mag = sqrt(a .^ 2 +  b .^ 2);
    
        ang = atan(b ./ a);
    else
        Phi = [cos(2 * pi * f .* t), sin(2 * pi * f .* t), ones(height(t), width(f))];

        Phi = Phi(:, reshape(reshape(1 : width(Phi), width(f), 3)', 1, width(Phi)));
    
        Phi = mat2cell(Phi, height(t), 3 * ones(1, width(f)));

        mag = zeros(width(f), 1);
        
        ang = zeros(width(f), 1);

        for i = 1 : width(f)
            res = Phi{i} \ X(:, i);
            
            mag(i) = sqrt(res(1) .^ 2 + res(2) .^ 2);

            ang(i) = atan(res(2) ./ res(1));
        end
    end