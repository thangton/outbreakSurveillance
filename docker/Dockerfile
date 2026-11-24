# Utility image: python helpers in bin/ + database downloaders (curl, NCBI datasets)
FROM python:3.11-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl unzip ca-certificates procps \
 && rm -rf /var/lib/apt/lists/* \
 && ARCH=$(dpkg --print-architecture | sed 's/amd64/linux-amd64/;s/arm64/linux-arm64/') \
 && curl -fsSL -o /usr/local/bin/datasets \
      https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/${ARCH}/datasets \
 && chmod +x /usr/local/bin/datasets

RUN pip install --no-cache-dir biopython==1.84 matplotlib==3.9.2

LABEL org.opencontainers.image.source="https://github.com/thangton/outbreakSurveillance"
