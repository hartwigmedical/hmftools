package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleRegionZipper.updateRegionsWithCopyNumbers;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.copynumber.freec.FreecFileLoader;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatioFactory;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatioRegions;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFactory;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFile;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.FittedRegionWriter;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionFactory;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.purple.segment.PurpleSegmentFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfSlicerFileLoader;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFGermlineFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.structural.StructuralVariantFileLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurityPloidyEstimateApplication {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    static final double MIN_PURITY_DEFAULT = 0.05;
    static final double MAX_PURITY_DEFAULT = 1.0;
    static final double MIN_NORM_FACTOR = 0.33;
    static final double MAX_NORM_FACTOR = 2.0;

    private static final double MIN_REF_ALLELE_FREQUENCY = 0.4;
    private static final double MAX_REF_ALLELE_FREQUENCY = 0.65;
    private static final int MIN_COMBINED_DEPTH = 10;
    private static final int MAX_COMBINED_DEPTH = 100;
    private static final int MAX_PLOIDY = 20;
    private static final double PURITY_INCREMENTS = 0.01;
    private static final double NORM_FACTOR_INCREMENTS = 0.01;

    private static final String MIN_PURITY = "min_purity";
    private static final String MAX_PURITY = "max_purity";
    private static final String DB_ENABLED = "db_enabled";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String VCF_EXTENSION = "vcf_extension";
    private static final String VCF_EXTENSION_DEFAULT = ".annotated_sliced.vcf";
    private static final String CNV_RATIO_WEIGHT_FACTOR = "cnv_ratio_weight_factor";
    private static final double CNV_RATIO_WEIGHT_FACTOR_DEFAULT = 0.2;

    private static final String STRUCTURAL_VCF_EXTENTION = "structural_variant_extension";
    private static final String STRUCTURAL_VCF_EXTENTION_DEFAULT = "somaticSV.vcf.gz";

    private static final String PLOIDY_PENALTY_EXPERIMENT = "ploidy_penalty_experiment";

    private static final String OBSERVED_BAF_EXPONENT = "observed_baf_exponent";
    private static final double OBSERVED_BAF_EXPONENT_DEFAULT = 1;

    public static void main(final String... args) throws ParseException, IOException, HartwigException, SQLException {
        new PurityPloidyEstimateApplication(args);
    }

    private PurityPloidyEstimateApplication(final String... args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        if (runDirectory == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Purity Ploidy Estimator (PURPLE)", options);
            System.exit(1);
        }

        LOGGER.info("Loading germline variant data");
        final String vcfExtension = defaultValue(cmd, VCF_EXTENSION, VCF_EXTENSION_DEFAULT);
        final VCFGermlineFile vcfFile = VCFFileLoader.loadGermlineVCF(runDirectory, vcfExtension);
        final List<GermlineVariant> variants = VariantFilter.passOnly(vcfFile.variants());
        final String refSample = vcfFile.refSample();
        final String tumorSample = vcfFile.tumorSample();

        LOGGER.info("Loading {} Freec data", tumorSample);
        final String freecDirectory = freecDirectory(cmd, runDirectory, refSample, tumorSample);
        final List<FreecRatio> tumorRatio = FreecRatioFactory.loadTumorRatios(freecDirectory, tumorSample);
        // KODU: Even though this retrieves normal ratios, freec uses the tumor sample name in the file name.
        final List<FreecRatio> normalRatio = FreecRatioFactory.loadNormalRatios(freecDirectory, tumorSample);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(tumorRatio);

        final String structuralVariantExtension = defaultValue(cmd, STRUCTURAL_VCF_EXTENTION, STRUCTURAL_VCF_EXTENTION_DEFAULT);
        final String structuralVariantFile = PathExtensionFinder.build().findPath(runDirectory, structuralVariantExtension).toString();
        LOGGER.info("Loading structural variants from {}", structuralVariantFile);
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantFile);

        LOGGER.info("Merging structural variants into freec segmentation");
        final List<PurpleSegment> segments = PurpleSegmentFactory.createSegments(regions, structuralVariants);

        LOGGER.info("Mapping all observations to the regions defined by the tumor ratios");
        final ObservedRegionFactory observedRegionFactory =
                new ObservedRegionFactory(MIN_REF_ALLELE_FREQUENCY, MAX_REF_ALLELE_FREQUENCY, MIN_COMBINED_DEPTH, MAX_COMBINED_DEPTH);
        final List<ObservedRegion> observedRegions = observedRegionFactory.combine(segments, variants, tumorRatio, normalRatio);

        final Gender gender = Gender.fromObservedRegions(observedRegions);
        LOGGER.info("Sample gender is {}", gender.toString().toLowerCase());
        final double cnvRatioWeight = defaultValue(cmd, CNV_RATIO_WEIGHT_FACTOR, CNV_RATIO_WEIGHT_FACTOR_DEFAULT);
        final boolean ploidyPenaltyExperiment = cmd.hasOption(PLOIDY_PENALTY_EXPERIMENT);
        final double observedBafExponent = defaultValue(cmd, OBSERVED_BAF_EXPONENT, OBSERVED_BAF_EXPONENT_DEFAULT);
        final FittedRegionFactory fittedRegionFactory =
                new FittedRegionFactory(gender, MAX_PLOIDY, cnvRatioWeight, ploidyPenaltyExperiment, observedBafExponent);

        LOGGER.info("Fitting purity");
        final double minPurity = defaultValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
        final double maxPurity = defaultValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(MAX_PLOIDY,
                minPurity,
                maxPurity,
                PURITY_INCREMENTS,
                MIN_NORM_FACTOR,
                MAX_NORM_FACTOR,
                NORM_FACTOR_INCREMENTS,
                fittedRegionFactory,
                observedRegions);

        Optional<FittedPurity> optionalBestFit = fittedPurityFactory.bestFit();
        if (optionalBestFit.isPresent()) {
            final FittedPurity bestFit = optionalBestFit.get();
            final List<FittedRegion> fittedRegions = fittedRegionFactory.fitRegion(bestFit.purity(), bestFit.normFactor(), observedRegions);

            final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, bestFit.purity(), bestFit.normFactor());
            final PurpleCopyNumberFactory purpleCopyNumberFactory = new PurpleCopyNumberFactory(purityAdjuster, fittedRegions);
            final List<PurpleCopyNumber> highConfidence = purpleCopyNumberFactory.highConfidenceRegions();
            final List<PurpleCopyNumber> smoothRegions = purpleCopyNumberFactory.smoothedRegions();

            final FittedPurityScore score = FittedPurityScoreFactory.score(fittedPurityFactory.allFits(), smoothRegions);
            final List<FittedRegion> enrichedFittedRegions = updateRegionsWithCopyNumbers(fittedRegions, highConfidence, smoothRegions);

            if (cmd.hasOption(DB_ENABLED)) {
                LOGGER.info("Persisting to database");
                final DatabaseAccess dbAccess = databaseAccess(cmd);
                dbAccess.writePurity(tumorSample, score, fittedPurityFactory.bestFitPerPurity());
                dbAccess.writeCopynumbers(tumorSample, smoothRegions);
                dbAccess.writeCopynumberRegions(tumorSample, enrichedFittedRegions);
                dbAccess.writeGeneCopynumberRegions(tumorSample, geneCopyNumbers(smoothRegions));
                dbAccess.writeStructuralVariants(tumorSample, structuralVariants);
            }

            final String outputDirectory = defaultValue(cmd, OUTPUT_DIRECTORY, runDirectory + File.separator + OUTPUT_DIRECTORY_DEFAULT);
            LOGGER.info("Writing to file location: {}", outputDirectory);

            PurpleCopyNumberFile.write(outputDirectory, tumorSample, smoothRegions);
            FittedPurityFile.write(outputDirectory, tumorSample, fittedPurityFactory.bestFitPerPurity());
            FittedPurityScoreFile.write(outputDirectory, tumorSample, score);
            FittedRegionWriter.writeCopyNumber(outputDirectory, tumorSample, enrichedFittedRegions);
        }

        LOGGER.info("Complete");
    }

    private List<GeneCopyNumber> geneCopyNumbers(List<PurpleCopyNumber> copyNumbers) throws IOException, EmptyFileException {
        final InputStream inputStream = getClass().getResourceAsStream("/bed/hmf_gene_panel.tsv");
        final List<HmfGenomeRegion> hmgGenomeRegions = HmfSlicerFileLoader.fromInputStream(inputStream);
        return GeneCopyNumberFactory.geneCopyNumbers(hmgGenomeRegions, copyNumbers);
    }

    @NotNull
    private static String defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, @NotNull final String defaultValue) {
        return cmd.hasOption(opt) ? cmd.getOptionValue(opt) : defaultValue;
    }

    private static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    @NotNull
    private static String freecDirectory(@NotNull final CommandLine cmd, @NotNull final String runDirectory,
            @NotNull final String refSample, @NotNull final String tumorSample) {
        return cmd.hasOption(FREEC_DIRECTORY)
                ? cmd.getOptionValue(FREEC_DIRECTORY)
                : FreecFileLoader.getFreecBasePath(runDirectory, refSample, tumorSample);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(OBSERVED_BAF_EXPONENT, true, "Observed baf exponent. Default 1");
        options.addOption(PLOIDY_PENALTY_EXPERIMENT, false, "Use experimental ploidy penality.");
        options.addOption(OUTPUT_DIRECTORY, true, "The output path. Defaults to freec_dir.");
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run.");
        options.addOption(FREEC_DIRECTORY, true, "The freec data path. Defaults to ../copyNumber/sampleR_sampleT/freec/");
        options.addOption(VCF_EXTENSION, true, "VCF file extension. Defaults to " + VCF_EXTENSION_DEFAULT);
        options.addOption(CNV_RATIO_WEIGHT_FACTOR, true, "CNV ratio deviation scaling.");

        options.addOption(MIN_PURITY, true, "Minimum purity (default 0.05)");
        options.addOption(MAX_PURITY, true, "Maximum purity (default 1.0)");

        options.addOption(DB_ENABLED, false, "Persist data to DB.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
