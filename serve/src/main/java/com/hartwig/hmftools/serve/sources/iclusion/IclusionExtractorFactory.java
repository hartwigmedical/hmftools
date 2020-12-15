package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;

public final class IclusionExtractorFactory {

    private IclusionExtractorFactory() {
    }

    @NotNull
    public static IclusionExtractor buildIclusionExtractor(@NotNull EventClassifierConfig config, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap,
            @NotNull DoidLookup missingDoidLookup) {
        return new IclusionExtractor(EventExtractorFactory.create(config, proteinResolver, driverGenes, allGenesMap),
                new ActionableTrialFactory(missingDoidLookup));
    }
}
