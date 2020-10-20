package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.curation.FusionCuration;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelEventExtractor.class);

    //TODO
//     "SPLICE VARIANT 7",
//             "Splice",
//             "DNMT3B7",
//             "LCS6-variant",
//             "AR-V7",
//             "ARv567es");

    private static final Set<String> DETAILLED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO = Sets.newHashSet("MUTATION",
            "mutant",
            "mut",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION", "FRAMESHIFT MUTATION");
    private static final Set<String> DETAILLED_GENE_LEVEL_INFO_WITH_TSG = Sets.newHashSet("inact mut",
            "biallelic inactivation",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "DELETERIOUS MUTATION",
            "negative");
    private static final Set<String> DETAILLED_GENE_LEVEL_INFO_WITH_ONCO = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "act mut",
            "ACTIVATING MUTATION",
            "Oncogenic Mutations",
            "pos",
            "positive");

    public GeneLevelEventExtractor() {
    }

    @VisibleForTesting
    @NotNull
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        String eventDescription = feature.description(); // TODO extract only event without gene
        if (DETAILLED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription) || feature.provenanceRule() != null) {
            for (DriverGene driverGene : driverGenes) {
                if (driverGene.gene().equals(feature.geneSymbol())) {
                    if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals("gene_only")) {
                                return GeneLevelEvent.ACTIVATION;
                            }
                        } else if (DETAILLED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription)) {
                            return GeneLevelEvent.ACTIVATION;
                        }
                    } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals("gene_only")) {
                                return GeneLevelEvent.INACTIVATION;
                            }
                        } else if (DETAILLED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription)) {
                            return GeneLevelEvent.ACTIVATION;
                        }
                    }
                }
            }
            LOGGER.warn("Gene {} is not present in driver catalog", feature.geneSymbol());
        } else if (DETAILLED_GENE_LEVEL_INFO_WITH_TSG.contains(eventDescription)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (DETAILLED_GENE_LEVEL_INFO_WITH_ONCO.contains(eventDescription)) {
            return GeneLevelEvent.ACTIVATION;
        }

        return GeneLevelEvent.UNKONWN;
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry,
            @NotNull List<DriverGene> driverGenes) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.GENE_LEVEL) {
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder()
                                .gene(feature.geneSymbol())
                                .event(extractGeneLevelEvent(feature, driverGenes))
                                .build());

            } else if (feature.type() == FeatureType.FUSION_PROMISCUOUS) {

                String curatedPromiscuousFusion = FusionCuration.curatedFusions(feature.geneSymbol());
                //TODO: check if this is needed
                // if (function.equals("Likely Loss-of-function")) {
                //            gene = Strings.EMPTY;
                //            typeEvent = Strings.EMPTY;
                //        }
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder().gene(curatedPromiscuousFusion).event(GeneLevelEvent.FUSION).build());
            }

        }
        return geneLevelEventsPerFeature;
    }

}
