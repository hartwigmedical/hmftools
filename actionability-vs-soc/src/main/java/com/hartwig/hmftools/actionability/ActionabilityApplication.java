package com.hartwig.hmftools.actionability;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.actionability.CNVs.ActionabilityCNVsAnalyzer;
import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeAnalyzer;
import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeMappingReading;
import com.hartwig.hmftools.actionability.variants.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityApplication {

    private static final Logger LOGGER = LogManager.getLogger(ActionabilityApplication.class);
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";
    private static final String PURPLE_DIRECTORY = "purple";

    public static void main(final String... args) throws ParseException, IOException {
        LOGGER.info("Determining actionability variants");
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);

        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDir);
        final List<SomaticVariant> variants = loadPassedSomaticVariants(run.tumorSample(), runDir);
        final List<GeneCopyNumber> geneCopyNumbers = loadPurpleGeneCopyNumbers(runDir, run.tumorSample());

        final String patientIdentifier = toPatientIdentifier(run.tumorSample());
        LOGGER.info("Tumor sample: " + run.tumorSample());
        LOGGER.info("patientId: " + patientIdentifier);

        // TODO: LISC
        // compare primary tumor locations with doid
        LOGGER.info("DOIDs of tumor location");
        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(cmd.getOptionValue(TUMOR_LOCATION_CSV));
        final PatientTumorLocation patientTumorLocation = extractPatientTumorLocation(patientTumorLocations, run.tumorSample());
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        final String cancerSubtype = patientTumorLocation != null ? patientTumorLocation.cancerSubtype() : Strings.EMPTY;
        LOGGER.info("Tumor location from patient: " + primaryTumorLocation);
        LOGGER.info("Cancer subtype from patient: " + cancerSubtype);

        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doids = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);
        LOGGER.info("DOID: " + doids);

        String fileCancerTumorsWithDOID = "/data/common/dbs/knowledgebases/output/knowledgebaseCancerTypes.tsv";

        LOGGER.info("");
        LOGGER.info("Start processing actionability somaticVariants");

        String fileActionabilityVariants = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
        String fileActionabilityRanges = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";

        LOGGER.info("Variants: " + variants.size());
        if (Files.exists(new File(fileActionabilityVariants).toPath()) && Files.exists(new File(fileActionabilityRanges).toPath())
                && Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
            ActionabilityVariantsAnalyzer analyzer =
                    ActionabilityVariantsAnalyzer.loadFromFileVariantsAndFileRanges(fileActionabilityVariants, fileActionabilityRanges);
            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(fileCancerTumorsWithDOID);

            Set<String> actionableGenes = analyzer.actionableGenes();
            LOGGER.info(actionableGenes.size() + " actionable genes found for variants (and ranges)");
            LOGGER.info("");

            List<SomaticVariant> variantsOnActionableGenes =
                    variants.stream().filter(variant -> actionableGenes.contains(variant.gene())).collect(Collectors.toList());

            LOGGER.info("Gene" + "\t" + "chromosome" + "\t" + "position" + "\t" + "ref" + "\t" + "alt" + "\t" + "drug" + "\t" + "drugsType"
                    + "\t" + "cancerType" + "\t" + "level" + "\t" + "response" + "\t" + "Actionable variant");
            for (SomaticVariant variant : variantsOnActionableGenes) {
                analyzer.actionableVariants(variant, cancerTypeAnalyzer, doids, primaryTumorLocation);
              //  analyzer.actionableRange(variant, cancerTypeAnalyzer, doids, primaryTumorLocation);
            }

        } else if (!Files.exists(new File(fileActionabilityVariants).toPath())) {
            LOGGER.warn("File does not exist: " + fileActionabilityVariants);
        } else if (!Files.exists(new File(fileActionabilityRanges).toPath())) {
            LOGGER.warn("File does not exist: " + fileActionabilityRanges);
        } else if (!Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
            LOGGER.warn("File does not exist: " + fileCancerTumorsWithDOID);
        }

        LOGGER.info("Finish processing actionability somaticVariants");
        LOGGER.info("");
        LOGGER.info("Start processing actionability cnvs");
        String fileActionabilityCNVs = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";

//        LOGGER.info("CNVs: " + geneCopyNumbers.size());
//        if (Files.exists(new File(fileActionabilityCNVs).toPath()) && Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            ActionabilityCNVsAnalyzer analyzerCNVs = ActionabilityCNVsAnalyzer.loadFromFileCNVs(fileActionabilityCNVs);
//            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(fileCancerTumorsWithDOID);
//            LOGGER.info("Gene" + "\t" + "cnvType" + "\t" + "drug" + "\t" + "drugType" + "\t" + "cancerType" + "\t"
//                    + "level" + "\t" + "response" + "\t" + "Actionable variant");
//            for (GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
//                analyzerCNVs.actionableCNVs(geneCopyNumber, cancerTypeAnalyzer, doids, primaryTumorLocation);
//            }
//        } else if (!Files.exists(new File(fileActionabilityCNVs).toPath())) {
//            LOGGER.warn("File does not exist: " + fileActionabilityCNVs);
//        } else if (!Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            LOGGER.warn("File does not exist: " + fileCancerTumorsWithDOID);
//        }
//        LOGGER.info("Finished processing actionability cnvs");
    }

    @NotNull
    private static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
        return GeneCopyNumberFile.read(fileName);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Complete path towards the (curated) tumor location csv.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
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
}
