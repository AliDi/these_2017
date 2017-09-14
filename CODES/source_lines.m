function source_lines(x_src ,F, nb_F)
    
    N=length(x_src);
    nb_F=length(F);
    for n=1:N
        hold on
        p1=plot(F,x_src(n)*ones(nb_F,1), 'color', [0.5 0.5 0.5],'linestyle','--');
        p1.Color(4) = 0.5;
    end