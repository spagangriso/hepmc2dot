test_cmd = py.test

install:
	pip install -r requirements.txt

test:
	$(test_cmd) --cov=.

testnocov:
	$(test_cmd)

.PHONY: install test testnocov
