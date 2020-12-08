package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.extraction.extractor.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.CopyNumberExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.ExonExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.FusionExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.GeneLevelExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.HotspotExtractor;
import com.hartwig.hmftools.serve.extraction.extractor.MutationTypeFilterAlgo;
import com.hartwig.hmftools.serve.extraction.extractor.SignatureExtractor;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccExtractorFactory {

    private static final Set<String> VALID_FUSION_GENES = Sets.newHashSet("IGH", "IGK", "IGL");

    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull EventClassifierConfig config, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap) {
        return buildViccExtractorWithInterpretationTsv(config, proteinResolver, driverGenes, allGenesMap, null);
    }

    @NotNull
    public static ViccExtractor buildViccExtractorWithInterpretationTsv(@NotNull EventClassifierConfig config,
            @NotNull ProteinResolver proteinResolver, @NotNull List<DriverGene> driverGenes,
            @NotNull Map<String, HmfTranscriptRegion> allGenesMap, @Nullable String featureInterpretationTsv) {
        Set<String> genesInExome = allGenesMap.keySet();
        GeneChecker exomeGeneChecker = new GeneChecker(genesInExome);

        Set<String> fusionGeneSet = Sets.newHashSet();
        fusionGeneSet.addAll(genesInExome);
        fusionGeneSet.addAll(VALID_FUSION_GENES);
        GeneChecker fusionGeneChecker = new GeneChecker(fusionGeneSet);

        MutationTypeFilterAlgo mutationTypeFilterAlgo = new MutationTypeFilterAlgo(driverGenes);
        return new ViccExtractor(new HotspotExtractor(exomeGeneChecker, proteinResolver, config.proteinAnnotationExtractor()),
                new CodonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, allGenesMap),
                new ExonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, allGenesMap),
                new GeneLevelExtractor(exomeGeneChecker,
                        fusionGeneChecker,
                        driverGenes,
                        config.activatingGeneLevelKeyPhrases(),
                        config.inactivatingGeneLevelKeyPhrases()),
                new CopyNumberExtractor(exomeGeneChecker),
                new FusionExtractor(fusionGeneChecker),
                new SignatureExtractor(config.microsatelliteUnstableEvents(),
                        config.highTumorMutationalLoadEvents(),
                        config.hrDeficiencyEvents()),
                featureInterpretationTsv);
    }
}
