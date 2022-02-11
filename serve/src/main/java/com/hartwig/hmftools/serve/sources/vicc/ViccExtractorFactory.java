package com.hartwig.hmftools.serve.sources.vicc;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;

import org.jetbrains.annotations.NotNull;

public final class ViccExtractorFactory {

    // For the current VICC release we have dealt with all driver inconsistencies we want to deal with, so want to ignore remaining.
    private static final String DEAL_WITH_DRIVER_INCONSISTENCIES_MODE = "filter";

    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull DoidLookup missingDoidLookup) {

        DealWithDriverInconsistentModeAnnotation annotation =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode(DEAL_WITH_DRIVER_INCONSISTENCIES_MODE);

        return new ViccExtractor(EventExtractorFactory.create(config, refGenomeResource, annotation),
                new ActionableEvidenceFactory(missingDoidLookup, new DrugCurator(), new EvidenceLevelCurator()));
    }
}
