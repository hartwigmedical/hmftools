package com.hartwig.hmftools.common.virus;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusInterpreterDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterDataLoader.class);

    private VirusInterpreterDataLoader() {
    }

    @NotNull
    public static VirusInterpreterData load(@NotNull String annotatedVirusTsv) throws IOException {
        LOGGER.info("Loading annotated virus data from {}", new File(annotatedVirusTsv).getParent());
        List<AnnotatedVirusV1> viruses = AnnotatedVirusFileV1.read(annotatedVirusTsv);

        List<AnnotatedVirusV1> reportable = Lists.newArrayList();
        List<AnnotatedVirusV1> unreported = Lists.newArrayList();
        for (AnnotatedVirusV1 virus : viruses) {
            if (virus.reported()) {
                reportable.add(virus);
            } else {
                unreported.add(virus);
            }
        }

        LOGGER.info(" Loaded {} annotated viruses (of which {} are reportable) from {}",
                viruses.size(),
                reportable.size(),
                annotatedVirusTsv);

        return ImmutableVirusInterpreterData.builder().unreportedViruses(unreported).reportableViruses(reportable).build();
    }
}
