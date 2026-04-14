FROM python:3.13

# 1. Force apt-get to NEVER ask for keyboard input (Fixes the infinite hang!)
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# 2. Install apt dependencies quietly
RUN apt-get update && apt-get install -y --no-install-recommends \
    mafft \
    clustalo \
    mrbayes \
    wget \
    tar \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 3. Download MUSCLE (Quietly to prevent log crashes)
RUN wget -q https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3 -O /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle

# 4. Download IQ-TREE (Removed '-v' from tar to prevent log viewer crash)
RUN wget -q https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-Linux.tar.gz \
    && tar -xzf iqtree-3.1.1-Linux.tar.gz \
    && mv iqtree-3.1.1-Linux/bin/* /usr/local/bin/ \
    && chmod +x /usr/local/bin/iqtree3 \
    && rm iqtree-3.1.1-Linux.tar.gz

# 5. Download MPBoot (Quiet download and extraction)
RUN wget -q https://github.com/diepthihoang/mpboot/releases/download/v1.2/mpboot-avx-1.2-Linux.tar.gz \
    && tar -xzf mpboot-avx-1.2-Linux.tar.gz \
    && mv mpboot-avx-1.2-Linux/mpboot-avx /usr/local/bin/mpboot \
    && chmod +x /usr/local/bin/mpboot \
    && rm mpboot-avx-1.2-Linux.tar.gz

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PORT=8080

CMD gunicorn main:app -b 0.0.0.0:$PORT -w 2 --timeout 300
