package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.codon.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.extraction.exon.ExonExtractor;
import com.hartwig.hmftools.serve.extraction.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.signature.SignatureExtractor;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;

import org.jetbrains.annotations.NotNull;

public final class EventExtractorFactory {

    private static final Set<String> VALID_FUSION_GENES = Sets.newHashSet("IGH", "IGK", "IGL");

    private EventExtractorFactory() {
    }

    @NotNull
    public static EventExtractor create(@NotNull EventClassifierConfig config, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap) {
        Set<String> genesInExome = allGenesMap.keySet();
        GeneChecker exomeGeneChecker = new GeneChecker(genesInExome);

        Set<String> fusionGeneSet = Sets.newHashSet();
        fusionGeneSet.addAll(genesInExome);
        fusionGeneSet.addAll(VALID_FUSION_GENES);
        GeneChecker fusionGeneChecker = new GeneChecker(fusionGeneSet);

        MutationTypeFilterAlgo mutationTypeFilterAlgo = new MutationTypeFilterAlgo(driverGenes);
        return new EventExtractor(new HotspotExtractor(exomeGeneChecker, proteinResolver, config.proteinAnnotationExtractor()),
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
                        config.hrDeficiencyEvents()));
    }
}
