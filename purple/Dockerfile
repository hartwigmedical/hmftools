FROM openjdk:17-alpine

RUN mkdir /usr/share/hartwig
ARG VERSION
ADD target/purple-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/purple.jar

ENTRYPOINT ["java", "-jar", "/usr/share/hartwig/purple.jar"]
