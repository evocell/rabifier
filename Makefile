.PHONY: clean test install

help:
	@echo "check-dependencies - check if the third party software is installed"
	@echo "test - run tests"
	@echo "clean - remove all temporary and build files"
	@echo "develop - install the package in the 'development mode'"
	@echo "install - install the package in the active Python's site-packages"
	@echo "build - build the package"
	@echo "publish - package and upload a release"
	@echo "publish-test - (test) package and upload a release"

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .eggs/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -exec rm -rf {} +

check-dependencies:
	@echo "Checking dependencies..."
	@for f in phmmer hmmbuild hmmpress hmmscan blastp meme mast cd-hit prank mafft; do \
		hash $$f 2>/dev/null && echo "  $$f (OK)" || echo "  $$f (Missing)"; \
	done

test: check-dependencies
	@echo "Running tests..."
	@python -m unittest discover --start-directory tests --verbose
	python setup.py test

develop: clean
	python setup.py develop

install: clean
	python setup.py install

build: clean
	python setup.py sdist bdist_wheel

publish: build 
	twine register -r pypitest dist/rabifier*.tar.gz
	twine upload -r pypi dist/*

publish-test: build 
	twine register -r pypitest dist/rabifier*.tar.gz
	twine upload -r pypitest dist/*
