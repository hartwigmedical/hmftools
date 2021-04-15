package com.hartwig.hmftools.serve.sources.iclusion;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class IclusionExtractorFactory {

    // For iClusion we want to deal with any driver inconsistency!
    private static final boolean REPORT_DRIVER_INCONSISTENCIES = true;

    private IclusionExtractorFactory() {
    }

    @NotNull
    public static IclusionExtractor buildIclusionExtractor(@NotNull EventClassifierConfig config,
            @NotNull RefGenomeResource refGenomeResource, @NotNull DoidLookup missingDoidLookup) {
        return new IclusionExtractor(EventExtractorFactory.create(config, refGenomeResource, REPORT_DRIVER_INCONSISTENCIES),
                new ActionableTrialFactory(missingDoidLookup));
    }
}
