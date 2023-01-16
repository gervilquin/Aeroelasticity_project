function [S, n] = compute_normals_surface(horse_shoe,c)
    S = zeros(length(horse_shoe.A),1);
    n = zeros(length(horse_shoe.A),3);
    for i=1:length(horse_shoe.A)
       S(i) = ((horse_shoe.C(i,2)-horse_shoe.B(i,2)))*c;
       %
       n_aux = -cross((horse_shoe.B(i,:)-horse_shoe.A(i,:)),...
           (horse_shoe.C(i,:)-horse_shoe.B(i,:)));
       n(i,:) = n_aux./norm(n_aux);
   end
end
