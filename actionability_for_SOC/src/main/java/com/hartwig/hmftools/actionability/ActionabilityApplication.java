package com.hartwig.hmftools.actionability;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.List;
import java.io.File;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.actionability.CNVs.ActionabilityCNVsAnalyzer;
import com.hartwig.hmftools.actionability.fusions.ActionabilityFusionAnalyzer;
import com.hartwig.hmftools.actionability.variants.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityApplication {


    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityApplication.class);
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";

    public static void main(final String... args) throws ParseException, IOException, SQLException {
        LOGGER.info("Determining actionability variants.");
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);

        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDir);
        final List<SomaticVariant> variants = loadPassedSomaticVariants(run.tumorSample(), runDir);
        final List<GeneCopyNumber> geneCopyNumbers;

        final String patientIdentifier = toPatientIdentifier(run.tumorSample());
        LOGGER.info("Tumor sample: " + run.tumorSample());
        LOGGER.info("patientId: " + patientIdentifier);

        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(cmd.getOptionValue(TUMOR_LOCATION_CSV));
        final PatientTumorLocation patientTumorLocation = extractPatientTumorLocation(patientTumorLocations, run.tumorSample());
        LOGGER.info("Tumor location from patient: " + patientTumorLocation.primaryTumorLocation());
        LOGGER.info("Cancer subtype from patient: " + patientTumorLocation.cancerSubtype());

        LOGGER.info("");
        LOGGER.info("Start processing actionability variants");

        String fileActionabilityVariants = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
        String fileActionabilityRanges = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";

        if (Files.exists(new File(fileActionabilityVariants).toPath()) && Files.exists(new File(fileActionabilityRanges).toPath())) {
            ActionabilityVariantsAnalyzer analyzer = ActionabilityVariantsAnalyzer.loadFromFileVariantsAndFileRanges(fileActionabilityVariants, fileActionabilityRanges);
            for (int i = 0; i < variants.size(); i ++) {
                LOGGER.info("Is actionable variant: " + analyzer.actionableVariants(variants.get(i), patientTumorLocation.primaryTumorLocation(), variants.size()));
                LOGGER.info("Is actionable ranges: " + analyzer.actionableRange(variants.get(i), patientTumorLocation.primaryTumorLocation(), variants.size()));

            }
        } else if (!Files.exists(new File(fileActionabilityVariants).toPath())){
            LOGGER.warn("File does not exist: " + fileActionabilityVariants);
        } else if(!Files.exists(new File(fileActionabilityRanges).toPath())){
            LOGGER.warn("File does not exist: " + fileActionabilityRanges);
        }

        LOGGER.info("");
        LOGGER.info("Start processing actionability cnvs");
        String fileActionabilityCNVs = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";

        if (Files.exists(new File(fileActionabilityCNVs).toPath())) {
            ActionabilityCNVsAnalyzer analyzerCNVs = ActionabilityCNVsAnalyzer.loadFromFileCNVs(fileActionabilityCNVs);
            for (int i = 0; i < variants.size(); i ++) {
                // change variants to gene copy number
                LOGGER.info("Is actionable CNVs: " + analyzerCNVs.actionableCNVs(variants.get(i), patientTumorLocation.primaryTumorLocation(), variants.size()));
            }
        } else {
            LOGGER.warn("File does not exist: " + fileActionabilityCNVs);
        }

        LOGGER.info("");
        LOGGER.info("Start processing actionability fusions");
        String fileActionabilityFusionPairs = "/data/common/dbs/knowledgebases/output/actionableFusionPairs.tsv";
        String fileActionabilityPromiscuousFive = "/data/common/dbs/knowledgebases/output/actionablePromiscuousFive.tsv";
        String fileActionabilityPromiscuousThree = "/data/common/dbs/knowledgebases/output/actionablePromiscuousThree.tsv";

        if (Files.exists(new File(fileActionabilityFusionPairs).toPath()) && Files.exists(new File(fileActionabilityPromiscuousFive).toPath())
        && Files.exists(new File(fileActionabilityPromiscuousThree).toPath())){
            ActionabilityFusionAnalyzer analyzerFusions = ActionabilityFusionAnalyzer.loadFromFileFusions(fileActionabilityFusionPairs,
                    fileActionabilityPromiscuousFive, fileActionabilityPromiscuousThree);
            for (int i = 0; i < variants.size(); i ++) {
                // change variants to structuralvariant, five_breakend, three_breakend, structuralvariantfusion
                LOGGER.info("Is actionable fusion: " + analyzerFusions.actionableFusions(variants.get(i), patientTumorLocation.primaryTumorLocation(), variants.size()));
            }
        } else if (!Files.exists(new File(fileActionabilityFusionPairs).toPath())){
            LOGGER.warn("File does not exist: " + fileActionabilityFusionPairs);
        } else if(!Files.exists(new File(fileActionabilityPromiscuousFive).toPath())){
            LOGGER.warn("File does not exist: " + fileActionabilityPromiscuousFive);
        } else if(!Files.exists(new File(fileActionabilityPromiscuousThree).toPath())){
            LOGGER.warn("File does not exist: " + fileActionabilityPromiscuousThree);
        }

        LOGGER.info("");
        LOGGER.info("Writing output data to file");

        LOGGER.info("");
        LOGGER.info("Finish processing actionability variants");
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
    private static List<SomaticVariant> loadPassedSomaticVariants(@NotNull final String sample, @NotNull final String path) throws IOException {
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
