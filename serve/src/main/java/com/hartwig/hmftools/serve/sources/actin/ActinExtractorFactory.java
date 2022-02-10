package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public class ActinExtractorFactory {

    private static final String DEAL_WITH_DRIVER_INCONSISTENCIES_MODE = "ignore";

    private ActinExtractorFactory() {
    }

    @NotNull
    public static ActinExtractor buildActinExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource) {
        DealWithDriverInconsistentModeAnnotation annotation =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode(DEAL_WITH_DRIVER_INCONSISTENCIES_MODE);
        return new ActinExtractor(EventExtractorFactory.create(config, refGenomeResource, annotation));
    }
}
