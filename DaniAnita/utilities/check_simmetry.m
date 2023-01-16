function check_simmetry(A,name)
    diff = max(max(abs(A-A')));
    if diff > 1e-2
        warning(strcat('Non symmetric<',name,'> diff = ',...
                    num2str(max(max(abs(A-A'))))));
    end
end