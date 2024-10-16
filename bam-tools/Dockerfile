FROM openjdk:17-alpine

ARG VERSION
ARG TOOL_NAME=bam-tools
ARG TOOL_ALIAS=bamtools
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

USER root
RUN apk update && apk upgrade && apk add \
    bash \
    procps

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/${TOOL_NAME}
RUN chmod +x /usr/local/bin/${TOOL_NAME}

RUN ln -s /usr/local/bin/${TOOL_NAME} /usr/local/bin/${TOOL_ALIAS}
