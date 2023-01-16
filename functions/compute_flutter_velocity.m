function flutter_vel = compute_flutter_velocity(p,u)

    flutter_vel = 0;
    for i = 1:length(u)-1
        if (p(i+1) >= 0) && (p(i) < 0)
            flutter_vel = -p(i)*(u(i+1)-u(i))/(p(i+1)-p(i)) + u(i);
        end
    end

end