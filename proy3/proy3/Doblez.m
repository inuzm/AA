function [ ps , caso ] = Doblez( B , g , delta)
% Funcion que obtiene una aproximaci?n al problmea de regi?n
% de confianza por el m?todo Doblez.


pN = -B \ g;
t = (g' * g) / (g' * B * g);
pC = -t * g;

p = -B\g;

if norm(p) <= delta
    
    caso = 'Newton';
    ps = p;
    
else
    pc = - ((g' * g) / (g' * B * g)) *g;
    if (norm(pc) >= delta)
        caso = 'Maximo Descenso';
        ps = -(delta / norm(g)) * g;
    else
        caso = 'Doblez';
        a = norm(p - pc)^2;
        b = 2 * ((p - pc)' * pc);
        c = pc' * pc - delta^2;
        disc = sqrt(b^2 - 4 * a * c);
        s1 = (-b + disc) / (2 * a);
        s2 = (-b - disc) / (2 * a);
        if 0 <= s1 <=1
            ps = pc + s1 * (p - pc);
        else
            ps = pc + s2 * (p - pc);
        end
    end
    
end
    
end

