FROM intel/oneapi-basekit:devel-ubuntu18.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y \
	&& apt-get install -y \
	git \
	gfortran \
	intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic \
	intel-oneapi-compiler-fortran

WORKDIR /tmp

# Install Atomistica - Python module and standalone code
RUN git clone https://github.com/Atomistica/atomistica.git \
	&& cd atomistica \
	&& cp setup.cfg.intel_omp setup.cfg \
	&& python3 setup.py install --prefix=/opt/atomistica --force \
	&& cd build_standalone \
	&& cp Makefile.intel Makefile \
	&& DEBUG=0 OPENMP=1 make factories \
	&& DEBUG=0 OPENMP=1 make mdcore \
	&& mkdir -p /opt/atomistica/bin \
	&& cp mdcore-* /opt/atomistica/bin


FROM ubuntu:18.04

RUN apt-get update -y \
	&& apt-get install -y \
	python3 \
	python3-pip \
	python3-numpy \
	python3-scipy

ENV PATH='/opt/atomistica/bin:/usr/sbin:/usr/bin:/sbin:/bin'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64'
ENV PYTHONPATH='/opt/atomistica/lib/python3.7/site-packages'

# Install Python modules
RUN python3 -m pip install ase

COPY --from=builder /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/ /opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin/
COPY --from=builder /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/ /opt/intel/oneapi/mkl/2021.3.0/lib/intel64/
COPY --from=builder /opt/atomistica/ /opt/atomistica/

CMD ["python3"]
