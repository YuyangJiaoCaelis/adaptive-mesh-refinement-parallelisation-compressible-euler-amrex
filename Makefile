AMREX_HOME ?= $(CURDIR)/../amrex

.PHONY: build check test clean

build:
	$(MAKE) -C solver/build -j8 AMREX_HOME="$(AMREX_HOME)"

check:
	find experiments -path '*/scripts/*.sh' -print0 | xargs -0 -n1 bash -n
	python3 -m compileall -q experiments/riemann_1d/analysis experiments/orientation_2d/analysis experiments/quadrant_scaling/analysis

test: check

clean:
	rm -rf solver/build/tmp_build_dir solver/build/main2d.gnu.MPI.ex
	find experiments -type d -name '__pycache__' -prune -exec rm -rf {} +
