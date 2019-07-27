function v = Rmag_to_e(v, delta_d, delta_i)
    
    R= [cos((pi*delta_d)/180)*cos((pi*delta_i)/180),	sin((pi*delta_d)/180),	cos((pi*delta_d)/180)*sin((pi*delta_i)/180);
		-sin((pi*delta_d)/180)*cos((pi*delta_i)/180),	cos((pi*delta_d)/180),	-sin((pi*delta_d)/180)*sin((pi*delta_i)/180);
		-sin((pi*delta_i)/180),	0,	cos((pi*delta_i)/180)];
    v = R*v;
    
end