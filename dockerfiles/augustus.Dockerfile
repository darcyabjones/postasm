ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

# Install htslib and samtools.
# Bam2wig in augustus depends on some intermediate samtools/htslib compilation
# rather than the actual headers/shared libraries, so I have to compile it
# separately.

RUN  apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       git \
       libbz2-dev \
       libcurl4-openssl-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_PREFIX="/opt/htslib"
ARG HTSLIB_VERSION="1.9"
ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"

ENV SAMTOOLS_PREFIX="/opt/samtools"
ARG SAMTOOLS_VERSION="1.9"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"

WORKDIR /tmp/htslib
RUN  git clone ${HTSLIB_REPO} . \
  && git fetch --tags \
  && git checkout ${HTSLIB_VERSION} \
  && autoheader \
  && autoconf \
  && ./configure --prefix="${HTSLIB_PREFIX}" --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo)

WORKDIR /tmp/samtools
RUN  git clone ${SAMTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout ${SAMTOOLS_VERSION} \
  && autoheader \
  && autoconf \
  && ./configure \
  && make -j $(grep -c ^processor /proc/cpuinfo)


# Install dependencies

RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       libbamtools-dev \
       libboost-all-dev \
       libboost-iostreams-dev \
       libboost-graph-dev \
       libcurl4-openssl-dev \
       libgsl-dev \
       liblpsolve55-dev \
       libssl-dev \
       libsuitesparse-dev \
       libsqlite3-dev \
       zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# This is for bam2wig
ENV TOOLDIR="/tmp"

ENV AUGUSTUS_PREFIX="/opt/augustus"
ARG AUGUSTUS_VERSION="8b1b14a7489e4545e89c8725dc33268f6c2a9117"
ARG AUGUSTUS_REPO="https://github.com/Gaius-Augustus/Augustus.git"


WORKDIR /tmp/augustus
RUN  git clone ${AUGUSTUS_REPO} . \
  && git fetch --tags \
  && git checkout ${AUGUSTUS_VERSION} \
  && mkdir bin \
  && sed -i "s/# SQLITE = true/SQLITE = true/g" common.mk \
  && sed -i "s/# COMPGENEPRED = true/COMPGENEPRED = true/g" common.mk \
  && sed -i 's~INSTALLDIR = .*~INSTALLDIR="${AUGUSTUS_PREFIX}"~g' Makefile \
  && cd auxprogs/bam2wig \
  && make \
  && cd /tmp/augustus \
  && make \
  && make install \
  && make test


FROM debian:${DEBIAN_VERSION}

ENV AUGUSTUS_PREFIX="/opt/augustus"
COPY --from=builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"

RUN  apt-get update \
  && apt-get install -y \
       libbamtools2.4.0 \
       libcurl3 \
       libgsl2 \
       libssl1.1 \
       libsqlite3-0 \
       lp-solve \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*


ENV PATH="${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

# This is useful for testing.
# COPY --from=builder "/tmp/augustus/examples" "${AUGUSTUS_PREFIX}/examples"
# RUN augustus --species=human --UTR=on ${AUGUSTUS_PREFIX}/examples/example.fa
