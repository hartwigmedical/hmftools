package com.hartwig.hmftools.common.vicc_knowledgebase;

import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VICC {

    private static final Logger LOGGER = LogManager.getLogger(VICC.class);

    @NotNull
    private final Map<String, FileBRCA> BRCA;

    VICC(@NotNull final Map<String, FileBRCA> BRCA) {
        this.BRCA = BRCA;
    }

}
