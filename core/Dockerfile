FROM ubuntu:latest

WORKDIR /fem

COPY . /fem/

RUN apt-get update && apt-get install -y gcc g++ make cmake

RUN apt-get install -y ligeigen3-dev libopenblas-dev liblapacke liblapacke-dev liblapacke-64 liblapacke-64-dev

RUN apt-get install -y gmsh libgmsh-dev

RUN apt-get install -y nlohmann-json3-dev

