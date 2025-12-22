FROM python:3.10-alpine
ARG VERSION=1.0.4

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1


#install gcc and other dependencies
#RUN apt-get update && \
#    apt-get install -y --no-install-recommends build-essential gcc && \
#    rm -rf /var/lib/apt/lists/*
#install python dependencies
#RUN pip install -U pip && \
#    pip install --no-cache-dir pandas==2.1.* && \
#    pip install --no-cache-dir scipy==1.11.* && \
#    pip install --no-cache-dir requests==2.31.*

#add trim kstar
COPY setup.py .
COPY kstar/ ./kstar/
RUN pip install --no-cache-dir . \
    && rm -rf ~/.cache/pip /var/lib/apt/lists/* 

#add resource files to /opt/kstar_resources/ and establish environment variable for resources
COPY RESOURCE_FILES/ /opt/kstar_resources/RESOURCE_FILES/
COPY NETWORKS/ /opt/kstar_resources/NETWORKS/
ENV KSTAR_RESOURCE_DIR=/opt/kstar_resources/

#install resource files
#COPY install_resources.py .
#RUN python install_resources.py

# Clean up pip cache and other unnecessary files
#RUN rm -rf ~/.cache/pip /var/lib/apt/lists/* 


#clean up unnecessary files
#RUN apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*