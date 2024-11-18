FROM bash:5.2

ARG VERSION
ARG TOOL_NAME=orange
ARG JAVA_DIR=/usr/share/java
ARG JAR=${JAVA_DIR}/${TOOL_NAME}_v${VERSION}.jar
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_NAME}

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}

## Required so that oncoanalyser can find the jar path. Oncoanalyser assumes:
## - the orange executable is a symlink to the run script
## - the jar is called orange.jar and is in the same dir as the run script
ADD target/entrypoint.sh ${JAVA_DIR}/entrypoint.sh
RUN ln -sf ${JAVA_DIR}/entrypoint.sh ${RUN_SCRIPT} && chmod +x ${RUN_SCRIPT}
RUN ln -sf ${JAR} ${JAVA_DIR}/${TOOL_NAME}.jar

RUN apk update && apk upgrade && apk add \
    procps \
    openjdk17

RUN ln -s /usr/local/bin/bash /bin/bash