FROM python:3.13
# help

RUN apt-get update && apt-get install -y wget tar \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
    mafft \
    clustalo \
    mrbayes \
    wget \
    tar \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# RUN wget https://github.com \
#     && chmod +x muscle-linux-x86.v5.3 \
#     && ./muscle-linux-x86.v5.3

# Download MUSCLE v5.3 (Linux x86_64)
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3 -O /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle

# Download IQ-TREE v3.1.1
RUN wget https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-Linux.tar.gz \
    && tar -xvzf iqtree-3.1.1-Linux.tar.gz \
    && mv iqtree-3.1.1-Linux/bin/* /usr/local/bin/ \
    && chmod +x /usr/local/bin/iqtree3

# Download MPBoot v1.2 (Corrected path)
RUN wget https://github.com/diepthihoang/mpboot/releases/download/v1.2/mpboot-avx-1.2-Linux.tar.gz \
    && tar -xvzf mpboot-avx-1.2-Linux.tar.gz \
    && mv mpboot-avx-1.2-Linux/mpboot-avx /usr/local/bin/mpboot \
    && chmod +x /usr/local/bin/mpboot

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PORT=8080

CMD gunicorn main:app -b 0.0.0.0:$PORT -w 2 --timeout 300
