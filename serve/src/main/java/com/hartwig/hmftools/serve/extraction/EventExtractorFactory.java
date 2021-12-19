package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicExtractor;
import com.hartwig.hmftools.serve.extraction.codon.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.extraction.exon.ExonExtractor;
import com.hartwig.hmftools.serve.extraction.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class EventExtractorFactory {

    private EventExtractorFactory() {
    }

    @NotNull
    public static EventExtractor create(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            boolean reportOnDriverInconsistencies) {
        Set<String> genesInExome = extractAllValidGenes(refGenomeResource.ensemblDataCache());
        GeneChecker exomeGeneChecker = new GeneChecker(genesInExome);

        Set<String> fusionGeneSet = Sets.newHashSet();
        fusionGeneSet.addAll(genesInExome);
        fusionGeneSet.addAll(extractAllGenesInvolvedInFusions(refGenomeResource.knownFusionCache()));
        GeneChecker fusionGeneChecker = new GeneChecker(fusionGeneSet);

        MutationTypeFilterAlgo mutationTypeFilterAlgo = new MutationTypeFilterAlgo(refGenomeResource.driverGenes());
        return new EventExtractor(new HotspotExtractor(exomeGeneChecker,
                refGenomeResource.proteinResolver(),
                config.proteinAnnotationExtractor()),
                new CodonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, refGenomeResource.ensemblDataCache()),
                new ExonExtractor(exomeGeneChecker, mutationTypeFilterAlgo, refGenomeResource.ensemblDataCache()),
                new GeneLevelExtractor(exomeGeneChecker,
                        fusionGeneChecker,
                        refGenomeResource.driverGenes(),
                        refGenomeResource.knownFusionCache(),
                        config.activatingGeneLevelKeyPhrases(),
                        config.inactivatingGeneLevelKeyPhrases(),
                        reportOnDriverInconsistencies),
                new CopyNumberExtractor(exomeGeneChecker, refGenomeResource.driverGenes(), reportOnDriverInconsistencies),
                new FusionExtractor(fusionGeneChecker,
                        refGenomeResource.knownFusionCache(),
                        config.exonicDelDupFusionKeyPhrases(),
                        reportOnDriverInconsistencies),
                new TumorCharacteristicExtractor(config.microsatelliteUnstableEvents(),
                        config.microsatelliteStableEvents(),
                        config.highTumorMutationalLoadEvents(),
                        config.lowTumorMutationalLoadEvents(),
                        config.hrDeficiencyEvents(),
                        config.hpvPositiveEvents(),
                        config.ebvPositiveEvents(),
                        config.hlaEvents()));
    }

    @NotNull
    private static Set<String> extractAllValidGenes(@NotNull EnsemblDataCache ensemblDataCache) {
        Set<String> genes = Sets.newHashSet();
        for (List<GeneData> genesPerChromosome : ensemblDataCache.getChrGeneDataMap().values()) {
            for (GeneData geneData : genesPerChromosome) {
                genes.add(geneData.GeneName);
            }
        }
        return genes;
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
