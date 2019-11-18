function plot_q2(q2)
	figure
	plot(real(q2),imag(q2));
	ax = gca;
	fb = ax.XBaseline;
	fb.Color = 'r';
	fb.BaseValue = 0;
	fb.Visible = 'on';
	fb = ax.YBaseline;
	fb.Color = 'r';
	fb.BaseValue = 0;
	fb.Visible = 'on';
	xlabel('Re(q^2)')
	ylabel('Im(q^2)')