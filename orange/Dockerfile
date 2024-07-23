FROM openjdk:17-alpine

RUN mkdir /usr/share/hartwig
ARG VERSION
ADD target/orange-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/orange.jar

#ENTRYPOINT ["java", "-jar", "/usr/share/hartwig/orange.jar"]
ENTRYPOINT ["find", "/"]
