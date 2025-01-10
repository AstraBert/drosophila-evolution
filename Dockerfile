FROM condaforge/miniforge3:latest

WORKDIR /app/

COPY conda_envs/*.yaml /app/conda_envs/

COPY setup.sh /app/

RUN bash /app/setup.sh

CMD ["/bin/bash"]

