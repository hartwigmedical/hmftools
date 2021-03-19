package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicExtractor;
import com.hartwig.hmftools.serve.extraction.codon.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.extraction.exon.ExonExtractor;
import com.hartwig.hmftools.serve.extraction.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;

import org.jetbrains.annotations.NotNull;

public final class EventExtractorFactory {

    private EventExtractorFactory() {
    }

    @NotNull
    public static EventExtractor create(@NotNull EventClassifierConfig config, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache,
            @NotNull Map<String, HmfTranscriptRegion> allGenesMap) {
        Set<String> genesInExome = allGenesMap.keySet();
        GeneChecker exomeGeneChecker = new GeneChecker(genesInExome);

        Set<String> fusionGeneSet = Sets.newHashSet();
        fusionGeneSet.addAll(genesInExome);
        fusionGeneSet.addAll(extractAllGenesInvolvedInFusions(knownFusionCache));
        GeneChecker fusionGeneChecker = new GeneChecker(fusionGeneSet);

        MutationTypeFilterAlgo mutationTypeFilterAlgo = new MutationTypeFilterAlgo(driverGenes);
        return new EventExtractor(new HotspotExtractor(exomeGeneChecker, proteinResolver, config.proteinAnnotationExtractor()),
                new CodonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, allGenesMap),
                new ExonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, allGenesMap),
                new GeneLevelExtractor(exomeGeneChecker,
                        fusionGeneChecker,
                        driverGenes,
                        knownFusionCache,
                        config.activatingGeneLevelKeyPhrases(),
                        config.inactivatingGeneLevelKeyPhrases()),
                new CopyNumberExtractor(exomeGeneChecker, driverGenes),
                new FusionExtractor(fusionGeneChecker, knownFusionCache),
                new TumorCharacteristicExtractor(config.microsatelliteUnstableEvents(),
                        config.highTumorMutationalLoadEvents(),
                        config.hrDeficiencyEvents(),
                        config.hpvPositiveEvents(),
                        config.ebvPositiveEvents()));
    }

    @NotNull
    private static Set<String> extractAllGenesInvolvedInFusions(@NotNull KnownFusionCache knownFusionCache) {
        Set<String> genes = Sets.newHashSet();
        for (KnownFusionData fusion : knownFusionCache.getData()) {
            if (fusion.Type == KnownFusionType.KNOWN_PAIR || fusion.Type == KnownFusionType.IG_KNOWN_PAIR
                    || fusion.Type == KnownFusionType.EXON_DEL_DUP) {
                genes.add(fusion.FiveGene);
                genes.add(fusion.ThreeGene);
            } else if (fusion.Type == KnownFusionType.PROMISCUOUS_5 || fusion.Type == KnownFusionType.IG_PROMISCUOUS) {
                genes.add(fusion.FiveGene);
            } else if (fusion.Type == KnownFusionType.PROMISCUOUS_3) {
                genes.add(fusion.ThreeGene);
            }
        }
        return genes;
    }
}
