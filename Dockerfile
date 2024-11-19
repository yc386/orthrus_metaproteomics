FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "orthrus", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
RUN echo "Make sure instanovo is installed:"
RUN python -c "import instanovo;print(instanovo.__version__)"

# The code to run when container is started:
COPY orthrus_v1 orthrus_v1
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "orthrus", "python", "orthrus_v1/annotated_orthrus_pt1.py"]