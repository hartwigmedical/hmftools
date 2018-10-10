package com.hartwig.hmftools.patientreporter.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariant;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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
//    private static final String FILE_ACTIONABILITY_CNVS = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static <T extends SomaticVariant> List<ActionabilityVariant> detectVariants(@NotNull List<T> variants,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;

        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        List<ActionabilityVariant> actionabilityVariants = Lists.newArrayList();
        if (Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath()) && Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())
                && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
            ActionabilityVariantsAnalyzer analyzer =
                    ActionabilityVariantsAnalyzer.loadFromFileVariantsAndFileRanges(FILE_ACTIONABILITY_VARIANTS, FILE_ACTIONABILITY_RANGES);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);

            Set<String> actionableGenesVariants = analyzer.actionableGenes();

            List<SomaticVariant> variantsOnActionableGenes =
                    variants.stream().filter(variant -> actionableGenesVariants.contains(variant.gene())).collect(Collectors.toList());

            for (SomaticVariant variant : variantsOnActionableGenes) {
                Set<ActionabilityVariant> data = analyzer.actionableVariants(variant, cancerTypeAnalyzer, doidsPrimaryTumorLocation);
                actionabilityVariants.addAll(data);
                LOGGER.info(data);
            }

        } else if (!Files.exists(new File(FILE_ACTIONABILITY_VARIANTS).toPath())) {
        LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_VARIANTS);
        } else if (!Files.exists(new File(FILE_ACTIONABILITY_RANGES).toPath())) {
        LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_RANGES);
        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
        LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
        }
        return actionabilityVariants;
    }
}
