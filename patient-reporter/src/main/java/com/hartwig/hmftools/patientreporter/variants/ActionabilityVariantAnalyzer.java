package com.hartwig.hmftools.patientreporter.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsEvidenceItems;
import com.hartwig.hmftools.common.actionability.fusion.ActionabilityFusionAnalyzer;
import com.hartwig.hmftools.common.actionability.fusion.ActionabilityFusionPairs;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceItems;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRangeEvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityVariantAnalyzer.class);
    private static final String FILE_CANCER_TUMORS_WITH_DOID = "/data/common/dbs/knowledgebases/output/knowledgebaseCancerTypes.tsv";
    private static final String FILE_ACTIONABILITY_VARIANTS = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
    private static final String FILE_ACTIONABILITY_RANGES = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";
    private static final String FILE_ACTIONABILITY_CNVS = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";
    private static final String FILE_ACTIONABILITY_FUSIONPAIRS = "/data/common/dbs/knowledgebases/output/actionableFusionPairs.tsv";
    private static final String FILE_ACTIONABILITY_PROMISCUOUS_FIVE =
            "/data/common/dbs/knowledgebases/output/actionablePromiscuousFive.tsv";
    private static final String FILE_ACTIONABILITY_PROMISCUOUS_THREE =
            "/data/common/dbs/knowledgebases/output/actionablePromiscuousThree.tsv";

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, VariantEvidenceItems> detectVariants(@NotNull List<T> variants,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        Map<T, VariantEvidenceItems> evidenceItemsPerVariant = Maps.newHashMap();
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        if (Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath()) && Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())
                && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            ActionabilityVariantsAnalyzer analyzer =
                    ActionabilityVariantsAnalyzer.loadFromFileVariantsAndFileRanges(FILE_ACTIONABILITY_VARIANTS, FILE_ACTIONABILITY_RANGES);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);

            Set<String> actionableGenesVariants = analyzer.actionableGenes();

            List<T> variantsOnActionableGenes =
                    variants.stream().filter(variant -> actionableGenesVariants.contains(variant.gene())).collect(Collectors.toList());

            for (T variant : variantsOnActionableGenes) {
                evidenceItemsPerVariant.put(variant, analyzer.actionableVariants(variant, cancerTypeAnalyzer, doidsPrimaryTumorLocation));
            }

        } else if (!Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_VARIANTS);
        } else if (!Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_RANGES);
        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
        }

        return evidenceItemsPerVariant;
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, ActionabilityRangeEvidenceItem> detectVariantsRanges(@NotNull List<T> variants,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        Map<T, ActionabilityRangeEvidenceItem> evidenceItemsPerVariantRange = Maps.newHashMap();
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        if (Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath()) && Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())
                && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            ActionabilityVariantsAnalyzer analyzer =
                    ActionabilityVariantsAnalyzer.loadFromFileVariantsAndFileRanges(FILE_ACTIONABILITY_VARIANTS, FILE_ACTIONABILITY_RANGES);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);

            Set<String> actionableGenesVariantsRanges = analyzer.actionableGenes();

            List<T> variantsOnActionableGenesRanges =
                    variants.stream().filter(variant -> actionableGenesVariantsRanges.contains(variant.gene())).collect(Collectors.toList());

            for (T variant : variantsOnActionableGenesRanges) {
                evidenceItemsPerVariantRange.put(variant, analyzer.actionableRange(variant, cancerTypeAnalyzer, doidsPrimaryTumorLocation));
            }

        } else if (!Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_VARIANTS);
        } else if (!Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_RANGES);
        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
        }

        return evidenceItemsPerVariantRange;
    }

    @NotNull
    public static <T extends GeneCopyNumber> Map<T, ActionabilityCNVsEvidenceItems> detectCNVs(@NotNull List<T> CNVs,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        Map<T, ActionabilityCNVsEvidenceItems> evidenceItemsPerVariantCNVs = Maps.newHashMap();
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        if (Files.exists(new File(FILE_ACTIONABILITY_CNVS).toPath()) && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            ActionabilityCNVsAnalyzer analyzerCNVs = ActionabilityCNVsAnalyzer.loadFromFileCNVs(FILE_ACTIONABILITY_CNVS);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);

            Set<String> actionableGenesCNVS = analyzerCNVs.actionableGenes();

            List<T> variantsOnActionableGenes = CNVs.stream()
                    .filter(geneCopyNumber -> actionableGenesCNVS.contains(geneCopyNumber.gene()))
                    .collect(Collectors.toList());

            for (T CNV : variantsOnActionableGenes) {
                evidenceItemsPerVariantCNVs.put(CNV, analyzerCNVs.actionableCNVs(CNV, cancerTypeAnalyzer, doidsPrimaryTumorLocation));
            }

        } else if (!Files.exists(new File(FILE_ACTIONABILITY_CNVS).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_CNVS);
        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
        }
        return evidenceItemsPerVariantCNVs;
    }

    @NotNull
    public static <T extends GeneFusion> Map<T, FusionEvidenceItems> detectFusions(@NotNull List<T> fusion,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        Map<T, FusionEvidenceItems> evidenceItemsPerVariantFusion = Maps.newHashMap();
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);
        if (Files.exists(new File(FILE_ACTIONABILITY_FUSIONPAIRS).toPath())
                && Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_FIVE).toPath()) && Files.exists(new File(
                FILE_ACTIONABILITY_PROMISCUOUS_THREE).toPath()) && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            ActionabilityFusionAnalyzer analyzerFusions = ActionabilityFusionAnalyzer.loadFromFileFusions(FILE_ACTIONABILITY_FUSIONPAIRS,
                    FILE_ACTIONABILITY_PROMISCUOUS_FIVE,
                    FILE_ACTIONABILITY_PROMISCUOUS_THREE);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);

//            Set<String> actionableGenesFusions = analyzerFusions.actionableGenes();
//            List<T> variantsOnActionableFusions = fusion.stream()
//                    .filter(fusionItem -> actionableGenesFusions.contains(fusionItem.reportable()))
//                    .collect(Collectors.toList());

//            for (T fusionGene : variantsOnActionableFusions) {
//                evidenceItemsPerVariantFusion.put(fusionGene, analyzerFusions.actionableFusions(cancerTypeAnalyzer, doidsPrimaryTumorLocation));
//            }


        } else if (!Files.exists(new File(FILE_ACTIONABILITY_FUSIONPAIRS).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_FUSIONPAIRS);
        } else if (!Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_FIVE).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_PROMISCUOUS_FIVE);
        } else if (!Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_THREE).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_PROMISCUOUS_THREE);
        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
        }
        return evidenceItemsPerVariantFusion;
    }
}
