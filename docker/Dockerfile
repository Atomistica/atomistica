FROM intel/oneapi:os-tools-ubuntu18.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y \
	&& apt-get install -y \
	curl \
	git \
	gfortran \
	python3 \
	python3-dev \
	python3-pip \
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic \
	intel-oneapi-compiler-fortran \
	intel-oneapi-mkl-devel

ENV PATH='/opt/intel/oneapi/mkl/2021.3.0/bin/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
ENV CPATH='/opt/intel/oneapi/mkl/2021.3.0/include:/opt/intel/oneapi/compiler/2021.3.0/linux/include'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/x64:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin'
ENV LIBRARY_PATH='/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2021.3.0/linux/lib'
ENV MKLROOT='/opt/intel/oneapi/mkl/2021.3.0'

# Numpy 1.19.5 and scipy 1.5.4 are last versions that support Ubuntu 18.04's Python 3.6
ENV NUMPY_VERSION="1.19.5"
ENV SCIPY_VERSION="1.5.4"

WORKDIR /tmp
RUN python3 -m pip install "setuptools<49.2.0" "cython==0.29.21" \
	&& curl -L https://github.com/numpy/numpy/archive/refs/tags/v${NUMPY_VERSION}.tar.gz | tar vxz
WORKDIR /tmp/numpy-${NUMPY_VERSION}
RUN python3 setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install

WORKDIR /tmp
RUN curl -L https://github.com/scipy/scipy/archive/refs/tags/v${SCIPY_VERSION}.tar.gz | tar vxz
WORKDIR /tmp/scipy-${SCIPY_VERSION}
RUN python3 setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install

WORKDIR /tmp
RUN git clone https://github.com/Atomistica/atomistica.git
WORKDIR /tmp/atomistica
# Install Atomistica - Python module and standalone code
RUN python3 -m pip install ase \
	&& cp setup.cfg.intel_omp setup.cfg \
	&& python3 setup.py install --prefix=/usr/local --force \
	&& cd build_standalone \
	&& cp Makefile.intel Makefile \
	&& DEBUG=0 OPENMP=1 make factories \
	&& DEBUG=0 OPENMP=1 make mdcore \
	&& mkdir -p /usr/local/bin \
	&& cp mdcore-* /usr/local/bin


FROM ubuntu:18.04

RUN apt-get update -y \
	&& apt-get install -y \
	python3 \
	python3-pip \
	strace \
	less

# Environment variables for MKL
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64'

# Copy numpy, scipy and atomistica
COPY --from=builder /usr/local/ /usr/local/

# Copy Fortran runtime - but only the libraries we actually need
COPY --from=builder /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/*.so* /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/

# Copy MKL - but only the libaries we actually need
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_rt.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_core.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_intel_thread.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_intel_lp64.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_def.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmkl_avx*.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/libmemkind.so* /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/

CMD ["python3"]
