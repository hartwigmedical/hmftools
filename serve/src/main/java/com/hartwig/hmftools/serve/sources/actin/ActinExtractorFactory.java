package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public class ActinExtractorFactory {

    // For ACTIN we don't care about driver inconsistencies
    private static final boolean REPORT_DRIVER_INCONSISTENCIES = false;

    private ActinExtractorFactory() {
    }

    @NotNull
    public static ActinExtractor buildActinExtractor(@NotNull EventClassifierConfig config,
            @NotNull RefGenomeResource refGenomeResource) {
        return new ActinExtractor(EventExtractorFactory.create(config, refGenomeResource, REPORT_DRIVER_INCONSISTENCIES));
    }
}
