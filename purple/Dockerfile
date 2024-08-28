FROM openjdk:17-slim

RUN mkdir /usr/share/hartwig
ARG VERSION
ADD target/purple-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/purple.jar
RUN apt update && apt install -y circos liblist-moreutils-perl

ENTRYPOINT ["find", "/resources"]
