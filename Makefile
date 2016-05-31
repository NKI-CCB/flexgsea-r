R=R --vanilla

all : test

test : test-r
.PHONY: test

test-r :
	echo 'requireNamespace("devtools"); devtools::test()' | $(R)
