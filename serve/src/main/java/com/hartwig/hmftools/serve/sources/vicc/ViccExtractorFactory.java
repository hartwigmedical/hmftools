package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccExtractorFactory {

    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull EventClassifierConfig config, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache,
            @NotNull Map<String, HmfTranscriptRegion> allGenesMap, @NotNull DoidLookup missingDoidLookup) {
        return buildViccExtractorWithInterpretationTsv(config,
                proteinResolver,
                driverGenes,
                knownFusionCache,
                allGenesMap,
                missingDoidLookup,
                null);
    }

    @NotNull
    public static ViccExtractor buildViccExtractorWithInterpretationTsv(@NotNull EventClassifierConfig config,
            @NotNull ProteinResolver proteinResolver, @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache,
            @NotNull Map<String, HmfTranscriptRegion> allGenesMap, @NotNull DoidLookup missingDoidLookup,
            @Nullable String featureInterpretationTsv) {
        return new ViccExtractor(EventExtractorFactory.create(config, proteinResolver, driverGenes, knownFusionCache, allGenesMap),
                new ActionableEvidenceFactory(missingDoidLookup),
                featureInterpretationTsv);
    }
}
