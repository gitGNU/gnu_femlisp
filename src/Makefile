all: TAGS

TAGS: *.lisp
	find . -name "*.lisp" -print | \
	etags - --regex='/[ \t]*\(definline\|defmemo\) +\([^ ]*\) //\1/'

linecount:
	find . -name "*.lisp" -exec wc -l '{}' ';'
	find . -name "*.lisp" -exec cat '{}' ';' | wc -l
	find . -name "*.lisp" -exec cat '{}' ';' | grep ";;" -v | wc -l
	find applications -name "*.lisp" -exec cat '{}' ';' | wc -l

clean:
	find . -name "*.x86f" -exec rm {} \;
	find . -name "*.err" -exec rm {} \;