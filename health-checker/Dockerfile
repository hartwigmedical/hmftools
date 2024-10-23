FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/health-checker
ADD target/health-checker-$VERSION-jar-with-dependencies.jar /opt/health-checker/health-checker.jar
