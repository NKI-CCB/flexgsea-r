R=R --vanilla

all : test
.PHONY: all

test : test-r check-r
.PHONY: test

test-r :
	echo 'requireNamespace("devtools"); devtools::test()' | $(R)
.PHONY: test-r
check-r :
	echo 'requireNamespace("devtools"); devtools::check()' | $(R)
.PHONY: check-r
fast-test-r :
	echo 'requireNamespace("devtools"); options(skip_expensive_tests=T); devtools::test()' | $(R)
.PHONY: fast-test-r
