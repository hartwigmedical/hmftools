FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/cider
ADD target/cider-$VERSION-jar-with-dependencies.jar /opt/cider/cider.jar
