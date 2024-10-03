# panasty_lacus_blast
Dockerized shiny blast server for _Astyanax lacustris_ pantranscriptome (Aciole Barbosa et a., 2024)

<br/>

Instructions to build the docker image from this git repo

note that a built image is already published at dockerhub:

https://hub.docker.com/r/davidaciole/panasty_lacus_blast

and a web app is available at shinyapps.io:

https://aciole-d.shinyapps.io/shiny_blast/

<br/>

### clone repo
```
git clone https://github.com/Aciole-David/panasty_lacus_blast.git
```
<br/>

### go to folder

```
cd panasty_lacus_blast/
```
<br/>

### decompress files within www/

```
gzip -d www/*.gz
```
<br/>

### build image

```
docker build -t panasty_lacus_blast .
```
<br/>

### run container exposing ports and using the built image
#note that it's possible to define the ports in the app.R and the Dockerfile

#but it must be the same in both files

#then you can change the '8080:8080' in the command below
```
docker run -p 8080:8080 panasty_lacus_blast
```
<br/>

### Navigate to the shiny interface page

open a web browser and go to either http://0.0.0.0:8080 or http://localhost:8080/









