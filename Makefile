all: www/index.html www/d3.v3.min.js

www/index.html: help.xhtml index.m4
	m4 -P index.m4 > www/index.html

help.xhtml: help.markdown
	pandoc -o help.xhtml help.markdown

www/d3.v3.min.js:
	wget http://d3js.org/d3.v3.min.js
