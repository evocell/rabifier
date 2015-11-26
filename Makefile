all:
	@:

test:
	@echo "Checking dependencies..."
	@echo "Running tests..."
	@python -m unittest discover --start-directory tests --verbose

