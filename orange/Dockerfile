FROM bash:5.2

ARG VERSION
ARG TOOL_NAME=orange
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_NAME}

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh ${RUN_SCRIPT}
RUN chmod +x ${RUN_SCRIPT}

RUN apk update && apk upgrade && apk add \
    procps \
    openjdk17

RUN ln -s /usr/local/bin/bash /bin/bash