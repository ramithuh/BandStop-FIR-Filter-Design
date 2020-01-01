
h_i=[ones(500,1); linspace(1,0,100).' ; zeros(1,300).' ; linspace(0,1,150).';ones(500,1)];

 
 
f = figure; 
p = plot(h_i(:,1));
p(1).LineWidth = 2;
