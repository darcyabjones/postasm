ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION}

ENV PYTHON3_PREFIX="/opt/python3"
ENV AUGUSTUS_PREFIX="/opt/augustus"

COPY --from="darcyabjones/augustus" "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from="darcyabjones/python3" "${PYTHON3_PREFIX}" "${PYTHON3_PREFIX}"


ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

ENV PATH="${PATH}:${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts:${PYTHON3_PREFIX}/bin"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${PYTHON3_PREFIX}/lib"
ENV LD_LIBRARY_PATH="${LIBRARY_PATH}"
ENV CPATH="${CPATH}:${PYTHON3_PREFIX}/include"

# ENV PYTHONHOME="/usr:/usr/local:${PYTHON3_PREFIX}:${PYTHONHOME}"
# ENV PYTHONPATH="/usr/bin:/usr/local/bin:${PYTHON3_PREFIX}/bin:${PYTHON3_PREFIX}/lib:${PYTHONPATH}"

ENV BUSCO_PREFIX="/opt/busco"
ENV BUSCO_REPO="https://gitlab.com/ezlab/busco.git"
ENV BUSCO_VERSION="1554283ab8ee7dd5b5290f4f748234f456c36e66"

RUN  apt-get update \
  && apt-get install -y \
       git \
       hmmer \
       libbamtools2.4.0 \
       libcurl3 \
       libgsl2 \
       libssl1.1 \
       libsqlite3-0 \
       lp-solve \
       ncbi-blast+ \
       perl \
       r-base \
       wget \
       zlib1g \
  && rm -rf /var/lib/apt/lists/* \
  && git clone "${BUSCO_REPO}" "${BUSCO_PREFIX}" \
  && cd "${BUSCO_PREFIX}" \
  && git fetch --tags \
  && git checkout "${BUSCO_VERSION}" \
  && python3 setup.py install \
  && sed -i 's~#!/usr/bin/env python~#!/usr/bin/env python3~g' scripts/run_BUSCO.py \
  && sed -i 's~#!/usr/bin/env python~#!/usr/bin/env python3~g' scripts/generate_plot.py \
  && cp config/config.ini.default config/config.ini \
  && sed -i "s~path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/bin/~path = ${AUGUSTUS_PREFIX}/bin/~g" config/config.ini \
  && sed -i "s~path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/scripts/~path = ${AUGUSTUS_PREFIX}/scripts/~g" config/config.ini \
  && sed -i 's~path = /home/osboxes/BUSCOVM/hmmer/hmmer-3.1b2-linux-intel-ia32/binaries/~path = /usr/bin/~g' config/config.ini

ENV PATH="${PATH}:${BUSCO_PREFIX}/scripts"
