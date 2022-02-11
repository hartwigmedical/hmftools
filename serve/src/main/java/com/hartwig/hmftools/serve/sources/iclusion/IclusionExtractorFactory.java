package com.hartwig.hmftools.serve.sources.iclusion;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class IclusionExtractorFactory {

    // For iClusion we want to deal with any driver inconsistency!
    private static final String DEAL_WITH_DRIVER_INCONSISTENCIES_MODE = "filter";

    private IclusionExtractorFactory() {
    }

    @NotNull
    public static IclusionExtractor buildIclusionExtractor(@NotNull EventClassifierConfig config,
            @NotNull RefGenomeResource refGenomeResource, @NotNull DoidLookup missingDoidLookup) {
        DealWithDriverInconsistentModeAnnotation annotation =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode(DEAL_WITH_DRIVER_INCONSISTENCIES_MODE);
        return new IclusionExtractor(EventExtractorFactory.create(config, refGenomeResource, annotation),
                new ActionableTrialFactory(missingDoidLookup));
    }
}
