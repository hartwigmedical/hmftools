package com.hartwig.hmftools.serve.sources.vicc.check;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneChecker {

    private static final Logger LOGGER = LogManager.getLogger(GeneChecker.class);

    private static final Set<String> GENES = Sets.newHashSet("IGH", "IGK", "IGL");

    public boolean isValidGene(@Nullable String gene, @Nullable HmfTranscriptRegion canonicalTranscript, @NotNull String event,
            @Nullable String eventExtractor) {

        if (canonicalTranscript != null) { //  Is gene part of canonical transcripts?
            return true;
        } else if (GENES.contains(gene)) { //  Is gene part of IG gene for fusions
            if (eventExtractor != null && eventExtractor.equals("fusion")) {
                return true;
            } else {
                LOGGER.warn("Could not find gene {} for event {} for fusion. Skipping extraction!", gene, event);
                return false;
            }
        } else {
            LOGGER.warn("Could not find gene {} for event {} as gene. Skipping extraction!", gene, event);
            return false;
        }
    }
}
