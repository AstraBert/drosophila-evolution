sudo docker ps -qf "name=drosophila_project_env1" > docker_check.log
if  [ -s "docker_check.log" ]
then
    echo "There is a container named 'drosophila_project_env1', deleting it to create a new one..."
    sudo docker rm --force $(sudo docker ps -qf "name=drosophila_project_env1")
else
    echo "There is no container named 'drosophila_project_env1', creating one..."
fi
rm -rf docker_check.log
sudo docker compose up -d
sudo docker exec -it $(sudo docker ps -qf "name=drosophila_project_env1") /bin/bash
