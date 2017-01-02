package com.hartwig.hmftools.patientreporter.slicing;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GenomeRegion {
    private static final Logger LOGGER = LogManager.getLogger(GenomeRegion.class);

    @NotNull
    private final String chromosome;
    private final long start;
    private final long end;

    public GenomeRegion(@NotNull final String chromosome, final long start, final long end) {
        if (end < start) {
            LOGGER.warn("Invalid region found: Chrom=" + chromosome + ", start=" + start + ", end=" + end);
        }
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;

    }

    long bases() {
        return 1 + end - start;
    }
}
