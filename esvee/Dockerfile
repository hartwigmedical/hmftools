FROM bash:5.2

ARG VERSION
ARG TOOL_NAME=esvee
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/${TOOL_NAME}
RUN chmod +x /usr/local/bin/${TOOL_NAME}

RUN apk update && apk upgrade && apk add \
    procps \
    openjdk17

RUN ln -s /usr/local/bin/bash /bin/bash
