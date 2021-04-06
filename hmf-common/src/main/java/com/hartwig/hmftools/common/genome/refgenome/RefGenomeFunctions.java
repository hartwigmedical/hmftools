package com.hartwig.hmftools.common.genome.refgenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class RefGenomeFunctions {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeFunctions.class);

    public static final String CHR_PREFIX = "chr";

    private RefGenomeFunctions() {
    }

    @NotNull
    public static String refGenomeChromosome(@NotNull final String chromosome, @NotNull RefGenomeVersion version) {
        if ((version == RefGenomeVersion.V38 || version == RefGenomeVersion.HG19) && !chromosome.contains(CHR_PREFIX)) {
            return CHR_PREFIX + chromosome;
        } else if (version == RefGenomeVersion.V37) {
            return stripChromosome(chromosome);
        } else {
            LOGGER.warn("Unrecognized ref genome version for making chromosome ref genome specific: {}", version);
            return chromosome;
        }
    }

    @NotNull
    public static String stripChromosome(@NotNull final String chromosome) {
        if (chromosome.startsWith(CHR_PREFIX)) {
            return chromosome.substring(CHR_PREFIX.length());
        }

        return chromosome;
    }
}
