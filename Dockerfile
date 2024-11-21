FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Install make and curl
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        make \
    && \
    apt-get update -y && \
    apt-get clean && \
    rm -rf /var/lib/{apt,dpkg,cache,log}


# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "orthrus", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
RUN echo "Make sure instanovo is installed:"
RUN python -c "import instanovo;print(instanovo.__version__)"

# The code to run when container is started:
COPY orthrus_v1 orthrus_v1
COPY Makefile Makefile
# ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "orthrus", "python", "orthrus_v1/annotated_orthrus_pt1.py"]
