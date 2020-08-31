FROM continuumio/miniconda3:4.8.2

LABEL description="Docker image for Deep Biosphere Geneflow project"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set temp dir
ENV TMPDIR="/tmp"

# Set workdir
WORKDIR /analysis

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && apt-get clean

# Add environment file
COPY environment.yaml .

# Install environment into base
RUN conda env update -n base -f environment.yaml && conda clean -a

# Set start up
CMD ["/bin/bash"]