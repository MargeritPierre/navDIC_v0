global hd
seed = hd.Seeds(2) ;
frame = hd.CurrentFrame ;

t = seed.Triangles ;
p = seed.MovingPoints(:,:,frame) ; seed.Points ;
f = seed.DataFields.u1(:,:,frame) ;
lvl = min(size(colormap,1),50)-1 ; linspace(0.2,1.6,10) ; logspace(log10(0.02),log10(2),20) ;

[C,lvl] = tricontours(t,p,f,lvl) ;

