package com.hartwig.hmftools.common.virus;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusInterpreterDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterDataLoader.class);

    @NotNull
    public static VirusInterpreterData load(@NotNull String annotatedVirusTsv) throws IOException
    {
        LOGGER.info("Loading VirusInterpreter data from {}", new File(annotatedVirusTsv).getParent());
        List<AnnotatedVirus> viruses = AnnotatedVirusFile.read(annotatedVirusTsv);

        LOGGER.info(" Loaded {} annotated viruses (of which {} are reportable) from {}",
                viruses.size(),
                viruses.stream().filter(AnnotatedVirus::reported).count(),
                annotatedVirusTsv);

        return ImmutableVirusInterpreterData.builder().allViruses(viruses).build();
    }
}
