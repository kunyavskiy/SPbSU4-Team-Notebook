beginfig(1)
	u := 1cm;
    path p;
    pair cur,v;
    p := (0,0);
    cur := (1u,0);
    v := (1u,0);
    for i := 1 upto 6:
    	p := p -- cur -- (0.5*u,sqrt(3)*u/2) -- cur;
    	v := v rotated 60;
    	cur := cur + v;
    endfor;
    p := p -- cycle;
    %draw p;
    for i := -50 upto 50:
    	for j := -50 upto 50:
    	    pair off;
    	    off := (i*1.5u,j*sqrt(3)*u+i*sqrt(3)*0.5u);
    		draw p shifted off withcolor white * 0.75;
    	endfor;
    endfor;	
endfig;

end.