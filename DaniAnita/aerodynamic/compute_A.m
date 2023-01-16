function A = compute_A(collocation_points,horse_shoe,n)
   for i=1:length(collocation_points)
        for j = 1:length(horse_shoe.A)
           vi12 = compute_induced_velocity(collocation_points(i,:),horse_shoe.A(j,:),horse_shoe.B(j,:));
           vi23 = compute_induced_velocity(collocation_points(i,:),horse_shoe.B(j,:),horse_shoe.C(j,:));
           vi34 = compute_induced_velocity(collocation_points(i,:),horse_shoe.C(j,:),horse_shoe.D(j,:));
           A(i,j) = dot(vi12+vi23+vi34,n(i,:));
        end
   end
end