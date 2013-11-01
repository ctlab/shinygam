www/index.html: www/help.xhtml www/index.m4
	(cd www; m4 -P index.m4 > index.html)

www/help.xhtml: www/help.markdown
	pandoc -o www/help.xhtml www/help.markdown
