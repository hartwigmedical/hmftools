package com.hartwig.hmftools.protect.virusinterpreter;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.protect.ProtectConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusInterpreterDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterDataLoader.class);

    private VirusInterpreterDataLoader() {
    }

    @NotNull
    public static List<AnnotatedVirus> load(@NotNull ProtectConfig config) throws IOException {
        return load(config.annotatedVirusTsv());
    }

    @NotNull
    public static List<AnnotatedVirus> load(@NotNull String annotatedVirusTsv) throws IOException {
        LOGGER.info("Loading annotated virus data from {}", new File(annotatedVirusTsv).getParent());
        List<AnnotatedVirus> annotatedViruses = AnnotatedVirusFile.read(annotatedVirusTsv);
        LOGGER.info(" Loaded {} annotated viruses from {}", annotatedViruses.size(), annotatedVirusTsv);
        return annotatedViruses;
    }
}
