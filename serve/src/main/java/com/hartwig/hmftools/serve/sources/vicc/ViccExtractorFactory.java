package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.serve.sources.vicc.extractor.CopyNumberExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.FusionExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneLevelExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneRangeExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.HotspotExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.SignatureExtractor;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccExtractorFactory {

    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull ProteinResolver proteinResolver, @NotNull List<DriverGene> driverGenes,
            @NotNull Map<String, HmfTranscriptRegion> allGenesMap) {
        return buildViccExtractorWithInterpretationTsv(proteinResolver, driverGenes, allGenesMap, null);
    }

    @NotNull
    public static ViccExtractor buildViccExtractorWithInterpretationTsv(@NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap,
            @Nullable String featureInterpretationTsv) {
        GeneChecker geneChecker = new GeneChecker(allGenesMap.keySet());

        return new ViccExtractor(new HotspotExtractor(proteinResolver, new ProteinAnnotationExtractor(), geneChecker),
                new CopyNumberExtractor(geneChecker),
                new FusionExtractor(geneChecker),
                new GeneLevelExtractor(geneChecker, driverGenes),
                new GeneRangeExtractor(geneChecker, driverGenes, allGenesMap),
                new SignatureExtractor(),
                featureInterpretationTsv);
    }
}
