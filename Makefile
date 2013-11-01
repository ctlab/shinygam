www/index.html: help.xhtml index.m4
	m4 -P index.m4 > www/index.html

help.xhtml: help.markdown
	pandoc -o help.xhtml help.markdown
