FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/patient-db
ADD target/patient-db-$VERSION-jar-with-dependencies.jar /opt/patient-db/patient-db.jar
