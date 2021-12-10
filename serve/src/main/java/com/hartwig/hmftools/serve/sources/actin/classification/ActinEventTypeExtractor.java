package com.hartwig.hmftools.serve.sources.actin.classification;

import java.util.List;

import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;
import com.hartwig.hmftools.serve.sources.actin.ActinTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ActinEventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ActinEventTypeExtractor.class);

    @NotNull
    private final GeneNameMapping37to38 geneNameMapping;

    public ActinEventTypeExtractor() {
        this.geneNameMapping = GeneNameMapping37to38.loadFromEmbeddedResource();
    }

    @NotNull
    public String extractGene(@NotNull List<ActinTrial> trials) {
        return Strings.EMPTY;
    }

    @NotNull
    public String extractEvent(@NotNull List<ActinTrial> trials) {
        return Strings.EMPTY;
    }
}
