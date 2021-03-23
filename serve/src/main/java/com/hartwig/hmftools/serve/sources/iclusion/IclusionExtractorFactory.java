package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class IclusionExtractorFactory {

    private IclusionExtractorFactory() {
    }

    @NotNull
    public static IclusionExtractor buildIclusionExtractor(@NotNull EventClassifierConfig config,
            @NotNull RefGenomeResource refGenomeResource, @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache,
            @NotNull DoidLookup missingDoidLookup) {
        return new IclusionExtractor(EventExtractorFactory.create(config, refGenomeResource, driverGenes, knownFusionCache),
                new ActionableTrialFactory(missingDoidLookup));
    }
}
