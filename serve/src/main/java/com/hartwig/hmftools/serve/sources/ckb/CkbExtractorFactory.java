package com.hartwig.hmftools.serve.sources.ckb;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class CkbExtractorFactory {

    // For CKB we want to ignore driver inconsistencies for now (their gene panel is larger than hmf driver panel to start with)
    private static final String DEAL_WITH_DRIVER_INCONSISTENCIES_MODE = "ignore";

    private CkbExtractorFactory() {
    }

    @NotNull
    public static CkbExtractor buildCkbExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource) {
        DealWithDriverInconsistentModeAnnotation annotation =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode(DEAL_WITH_DRIVER_INCONSISTENCIES_MODE);

        return new CkbExtractor(EventExtractorFactory.create(config, refGenomeResource, annotation),
                new ActionableEntryFactory());
    }
}
