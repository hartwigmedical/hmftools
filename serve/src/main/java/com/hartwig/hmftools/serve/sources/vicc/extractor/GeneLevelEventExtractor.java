package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
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

    public GeneLevelEventExtractor() {
    }

    @VisibleForTesting
    @NotNull
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        //                        "Oncogenic Mutations",
        //                        "TRUNCATING MUTATION",
        //                        "Truncating Mutations",
        //                        "DELETERIOUS MUTATION",
        //                        "feature_truncation",
        //                        "FRAMESHIFT TRUNCATION",
        //                        "FRAMESHIFT MUTATION",
        //                        "SPLICE VARIANT 7",
        //                        "Splice",
        //                        "DNMT3B7",
        //                        "LCS6-variant",
        //                        "AR-V7",
        //                        "ARv567es"
        if (feature.biomarkerType() != null) {
            if (feature.biomarkerType().equals("act mut")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("pos")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("positive")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("ACTIVATING MUTATION")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("Gain-of-function Mutations")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("Gain-of-Function")) {
                return GeneLevelEvent.ACTIVATION;
            } else if (feature.biomarkerType().equals("inact mut")) {
                return GeneLevelEvent.INACTIVATION;
            } else if (feature.biomarkerType().equals("biallelic inactivation")) {
                return GeneLevelEvent.INACTIVATION;
            } else if (feature.biomarkerType().equals("negative")) {
                return GeneLevelEvent.INACTIVATION;
            } else if (feature.biomarkerType().equals("Loss Of Function Variant")) {
                return GeneLevelEvent.INACTIVATION;
            } else if (feature.biomarkerType().equals("Loss Of Heterozygosity")) {
                return GeneLevelEvent.INACTIVATION;
            } else {
                for (DriverGene driverGene : driverGenes) {
                    if (driverGene.gene().equals(feature.geneSymbol())) {
                        if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                            if (feature.provenanceRule() != null) {
                                if (feature.provenanceRule().equals("gene_only")) {
                                    return GeneLevelEvent.ACTIVATION;
                                }
                            } else if (feature.biomarkerType() != null) {

                                if (feature.biomarkerType().equals("MUTATION")) {
                                    return GeneLevelEvent.ACTIVATION;
                                } else if (feature.biomarkerType().equals("mutant")) {
                                    return GeneLevelEvent.ACTIVATION;
                                } else if (feature.biomarkerType().equals("mut")) {
                                    return GeneLevelEvent.ACTIVATION;
                                }
                            }
                        } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                            if (feature.provenanceRule() != null) {
                                if (feature.provenanceRule().equals("gene_only")) {
                                    return GeneLevelEvent.INACTIVATION;
                                }
                            } else if (feature.biomarkerType() != null) {

                                if (feature.biomarkerType().equals("MUTATION")) {
                                    return GeneLevelEvent.ACTIVATION;
                                } else if (feature.biomarkerType().equals("mutant")) {
                                    return GeneLevelEvent.ACTIVATION;
                                } else if (feature.biomarkerType().equals("mut")) {
                                    return GeneLevelEvent.ACTIVATION;
                                }
                            }
                        }
                    }
                }
                LOGGER.warn("Gene {} is not present in driver catalog", feature.geneSymbol());
            }
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
