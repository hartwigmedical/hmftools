FROM openjdk:17-alpine

RUN mkdir /usr/share/hartwig
ARG VERSION
ADD target/virus-interpreter-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/virus-interpreter.jar

ENTRYPOINT ["java", "-jar", "/usr/share/hartwig/virus-interpreter.jar"]
