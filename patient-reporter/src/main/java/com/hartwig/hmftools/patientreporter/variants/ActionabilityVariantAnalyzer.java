package com.hartwig.hmftools.patientreporter.variants;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariant;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityVariantAnalyzer.class);
    private static String FILE_CANCER_TUMORS_WITH_DOID = "/data/common/dbs/knowledgebases/output/knowledgebaseCancerTypes.tsv";
    private static String FILE_ACTIONABILITY_VARIANTS = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
    private static String FILE_ACTIONABILITY_RANGES = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";
    private static String FILE_ACTIONABILITY_CNVS = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";
    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";
    private static final String PURPLE_DIRECTORY = "purple";


    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static List<ActionabilityVariant> detectVariants(@NotNull String runDirectory,
            @NotNull List<PatientTumorLocation> patientTumorLocations) throws IOException {

        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        final List<SomaticVariant> variants = loadPassedSomaticVariants(run.tumorSample(), runDirectory);

        final PatientTumorLocation patientTumorLocation = extractPatientTumorLocation(patientTumorLocations, runDirectory);
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

    @NotNull
    private static List<SomaticVariant> loadPassedSomaticVariants(@NotNull final String sample, @NotNull final String path)
            throws IOException {
        // TODO (KODU): Clean up once pipeline v3 no longer exists
        Path vcfPath;
        try {
            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V3);
        } catch (FileNotFoundException exception) {
            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V4);
        }
        return SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, vcfPath.toString());
    }

    @Nullable
    private static PatientTumorLocation extractPatientTumorLocation(@NotNull final List<PatientTumorLocation> patientTumorLocations,
            @NotNull final String sample) {
        final String patientIdentifier = toPatientIdentifier(sample);

        final List<PatientTumorLocation> matchingIdTumorLocations = patientTumorLocations.stream()
                .filter(patientTumorLocation -> patientTumorLocation.patientIdentifier().equals(patientIdentifier))
                .collect(Collectors.toList());

        // KODU: We should never have more than one curated tumor location for a single patient.
        assert matchingIdTumorLocations.size() < 2;

        if (matchingIdTumorLocations.size() == 1) {
            return matchingIdTumorLocations.get(0);
        } else {
            LOGGER.warn("Could not find patient " + patientIdentifier + " in clinical data!");
            return null;
        }
    }

    @NotNull
    private static String toPatientIdentifier(@NotNull final String sample) {
        if (sample.length() >= 12 && (sample.startsWith("CPCT") || sample.startsWith("DRUP"))) {
            return sample.substring(0, 12);
        }
        // KODU: If we want to generate a report for non-CPCT/non-DRUP we assume patient and sample are identical.
        return sample;
    }

    @NotNull
    private static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
        return GeneCopyNumberFile.read(fileName);
    }
}
