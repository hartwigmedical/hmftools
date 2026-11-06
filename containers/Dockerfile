# syntax=docker/dockerfile:1@sha256:87999aa3d42bdc6bea60565083ee17e86d1f3339802f543c0d03998580f9cb89

FROM ghcr.io/prefix-dev/pixi:0.73.0@sha256:66ae514a4421ca06587291beb037e5151b1b280a4f54eed82f57e9b70b94d5ca AS environment

ARG TOOL_NAME
ARG VERSION

WORKDIR /opt/hmftools
COPY pixi.toml pixi.lock ./

RUN --mount=type=cache,target=/root/.cache/rattler,sharing=locked \
    pixi install --locked --environment "${TOOL_NAME}"

# CUPPA's Python package is embedded in its JAR. Its dependencies are locked in
# pixi.lock; only the local package itself is installed here.
COPY --from=tool /target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar /tmp/tool.jar
RUN if [ "${TOOL_NAME}" = "cuppa" ]; then \
        pixi run --locked --environment cuppa unzip -q /tmp/tool.jar -d /tmp/cuppa \
        && pixi run --locked --environment cuppa python -m pip install --no-build-isolation --no-deps /tmp/cuppa/pycuppa; \
    fi \
    && rm /tmp/tool.jar \
    && printf '#!/bin/bash\n%s\nexport PATH="${JAVA_HOME}/bin:${PATH}"\nexec "$@"\n' "$(pixi shell-hook --environment "${TOOL_NAME}" --shell bash)" > /opt/hmftools/entrypoint \
    && chmod 0755 /opt/hmftools/entrypoint

FROM eclipse-temurin:17-jre-noble@sha256:543aebd60ff1deb9e906a8d4b117a7eda68a7f8e0d71041db2b5839d7fa057b8 AS runtime

ARG TOOL_NAME
ARG TOOL_ALIAS
ARG VERSION

WORKDIR /opt/hmftools
COPY --from=environment /opt/hmftools/.pixi/envs/${TOOL_NAME} /opt/hmftools/.pixi/envs/${TOOL_NAME}
COPY --from=environment /opt/hmftools/entrypoint /opt/hmftools/entrypoint
COPY --from=tool /target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar /usr/share/java/${TOOL_NAME}_v${VERSION}.jar
COPY --from=tool --chmod=0755 /target/entrypoint.sh /usr/share/java/run-tool

RUN ln -s /usr/share/java/run-tool "/usr/local/bin/${TOOL_ALIAS}" \
    && if [ "${TOOL_NAME}" = "orange" ]; then \
        ln -s "/usr/share/java/${TOOL_NAME}_v${VERSION}.jar" /usr/share/java/orange.jar; \
    fi

ENTRYPOINT ["/opt/hmftools/entrypoint"]
CMD ["bash"]
