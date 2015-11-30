help:
	@echo "check-dependencies - check if the third party software is installed"
	@echo "test - run tests"
	@echo "clean - remove all temporary and build files"
	@echo "install - install package"
	@echo "develop - install package in 'development mode'"
	@echo "publish - build and release the package"

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .eggs/
	find . -name "*.pyc" -exec rm -f {} +
	find . -name "__pycache__" -exec rm -rf {} +

check-dependencies:
	@echo "Checking dependencies..."
	@for f in phmmer hmmbuild hmmpress hmmscan blastp meme mast cd-hit prank mafft; do \
		hash $$f 2>/dev/null && echo "  $$f (OK)" || echo "  $$f (Missing)"; \
	done

test: check-dependencies
	@echo "Running tests..."
	@python -m unittest discover --start-directory tests --verbose

install: clean
	python setup.py install

develop: clean
	python setup.py develop

publish: clean
	python setup.py register
	python setup.py sdist upload
