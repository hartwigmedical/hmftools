package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
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

    private static final Set<String> VALID_FUSION_GENES = Sets.newHashSet("IGH", "IGK", "IGL");

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
        Set<String> genesInExome = allGenesMap.keySet();
        GeneChecker exomeGeneChecker = new GeneChecker(genesInExome);

        Set<String> fusionGeneSet = Sets.newHashSet();
        fusionGeneSet.addAll(genesInExome);
        fusionGeneSet.addAll(VALID_FUSION_GENES);
        GeneChecker fusionGeneChecker = new GeneChecker(fusionGeneSet);

        return new ViccExtractor(new HotspotExtractor(proteinResolver, new ProteinAnnotationExtractor(), exomeGeneChecker),
                new CopyNumberExtractor(exomeGeneChecker),
                new FusionExtractor(fusionGeneChecker),
                new GeneLevelExtractor(exomeGeneChecker, fusionGeneChecker, driverGenes),
                new GeneRangeExtractor(exomeGeneChecker, driverGenes, allGenesMap),
                new SignatureExtractor(),
                featureInterpretationTsv);
    }
}
