all: www/d3.v3.min.js

www/d3.v3.min.js:
	wget http://d3js.org/d3.v3.min.js -O www/d3.v3.min.js
