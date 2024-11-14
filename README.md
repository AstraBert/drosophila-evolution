# Unraveling the shared evolutionary history of European humans and Drosophila melanogaster
## Bachelor's thesis project by Astra Bertelli

**Supervisors**: Martin Kapun, Lino Ometto

**Collaborators**: Alan Bergland, Arnoud Estoup, Mathieu Gautier, Joaquin Nunez

This repository collects the data, the code and the analysis workflow for the project.

To reproduce the analysis, we strongly advise that you use the Docker image that we purposedly built for this project. You should get Docker for your platform following [this tutorial](https://dev.to/astrabert/1mindocker-2-get-docker-kh), adjust the path to the data you want to mount inside your Docker image in [`env`](./.env) under `ÙSERDATA_PATH` and then simply run:

```bash
docker compose up -d
docker exec -it $(sudo docker ps -qf "name=drosophila_project_env") /bin/bash
```

Or, if you are on Linux/macOS:

```bash
bash compose.sh
```

You'll be put inside the container (a semi-isolated virtual machine) in which our Docker image is up and running, and you'll find all your mounted data under `/gatk_modified/drosophila-project`.

If you want to check new releases for the image as well as past versions, you can do that on [Docker Hub](https://hub.docker.com/repository/docker/astrabert/silly-gat-kay/general) or on the [GitHub repository](https://github.com/AstraBert/silly-gat-kay)