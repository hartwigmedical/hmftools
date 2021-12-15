package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public class ActinExtractorFactory {

    // For iClusion we want to deal with any driver inconsistency!
    // For CKB we want to ignore driver inconsistencies for now (their gene panel is larger than hmf driver panel to start with)
    private static final boolean REPORT_DRIVER_INCONSISTENCIES = false;

    private ActinExtractorFactory() {
    }

    @NotNull
    public static ActinExtractor buildActinExtractor(@NotNull EventClassifierConfig config,
            @NotNull RefGenomeResource refGenomeResource) {
        return new ActinExtractor(EventExtractorFactory.create(config, refGenomeResource, REPORT_DRIVER_INCONSISTENCIES),
                new ActinTrialFactory());
    }
}
