FROM debian:jessie
MAINTAINER Jetstack <docker@jetstack.io>

ENV HOME /mykrobe

RUN adduser --disabled-password --gecos '' mykrobe && \
    mkdir $HOME && \
    chown -R mykrobe:mykrobe $HOME

COPY bin/* .$HOME/
COPY ./docker-entrypoint.sh $HOME/

USER mykrobe
WORKDIR $HOME

VOLUME ["$HOME/data"]

ENTRYPOINT ["./docker-entrypoint.sh"]
