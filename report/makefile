all:course_work.pdf

course_work.pdf: course_work.tex course_work.sty graphs/comp.tex graphs/fill.tex graphs/gnuplot-lua-tikz.sty
	latexmk -pdf course_work.tex
graphs/comp.tex:graphs/gen_comp.gp graphs/comp.data
	cd graphs; gnuplot gen_comp.gp

graphs/fill.tex:graphs/gen_fill.gp graphs/fill.data
	cd graphs; gnuplot gen_fill.gp
graphs/gnuplot-lua-tikz.sty:
	cd graphs; lua /usr/share/gnuplot/5.4/lua/gnuplot-tikz.lua style 
