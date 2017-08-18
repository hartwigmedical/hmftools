package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFactory.polyclonalProproption;
import static com.hartwig.hmftools.purple.PurpleRegionZipper.updateRegionsWithCopyNumbers;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.copynumber.freec.FreecGCContentFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFactory;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.BestFitFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.ratio.GCContent;
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
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.baf.BAFSupplier;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.SomaticConfig;
import com.hartwig.hmftools.purple.config.StructuralVariantConfig;
import com.hartwig.hmftools.purple.plot.ChartWriter;
import com.hartwig.hmftools.purple.ratio.ChromosomeLengthSupplier;
import com.hartwig.hmftools.purple.ratio.FreecRatioSupplier;
import com.hartwig.hmftools.purple.ratio.RatioSupplier;
import com.hartwig.hmftools.purple.ratio.ReadCountRatioSupplier;
import com.hartwig.hmftools.purple.segment.FreecSegmentSupplier;
import com.hartwig.hmftools.purple.segment.PCFSegmentSupplier;
import com.hartwig.hmftools.purple.somatic.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.purple.structural.StructuralVariantFileLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.Rserve.RserveException;

public class PurityPloidyEstimateApplication {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    static final double MIN_PURITY_DEFAULT = 0.08;
    static final double MAX_PURITY_DEFAULT = 1.0;
    static final double MIN_NORM_FACTOR_DEFAULT = 0.33;
    static final double MAX_NORM_FACTOR_DEFAULT = 2.0;

    private static final int MAX_PLOIDY = 20;
    private static final double PURITY_INCREMENTS = 0.01;
    private static final double NORM_FACTOR_INCREMENTS = 0.01;
    private static final int THREADS_DEFAULT = 2;

    private static final String COBALT = "cobalt";
    private static final String THREADS = "threads";
    private static final String MIN_PURITY = "min_purity";
    private static final String MAX_PURITY = "max_purity";
    private static final String MIN_NORM_FACTOR = "min_norm_factor";
    private static final String MAX_NORM_FACTOR = "max_norm_factor";
    private static final String DB_ENABLED = "db_enabled";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String GC_PROFILE = "gc_profile";
    private static final String PLOT = "plot";

    private static final String CNV_RATIO_WEIGHT_FACTOR = "cnv_ratio_weight_factor";
    private static final double CNV_RATIO_WEIGHT_FACTOR_DEFAULT = 0.2;

    private static final String PLOIDY_PENALTY_EXPERIMENT = "ploidy_penalty_experiment";

    private static final String OBSERVED_BAF_EXPONENT = "observed_baf_exponent";
    private static final double OBSERVED_BAF_EXPONENT_DEFAULT = 1;

    public static void main(final String... args)
            throws ParseException, IOException, HartwigException, SQLException, REXPMismatchException, RserveException, ExecutionException,
            InterruptedException {
        new PurityPloidyEstimateApplication(args);
    }

    private PurityPloidyEstimateApplication(final String... args)
            throws ParseException, IOException, HartwigException, SQLException, REXPMismatchException, RserveException, ExecutionException,
            InterruptedException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final int threads = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : THREADS_DEFAULT;
        final ExecutorService executorService = Executors.newFixedThreadPool(threads);
        try {
            // JOBA: Get common config
            final ConfigSupplier configSupplier = new ConfigSupplier(cmd, options);
            final CommonConfig config = configSupplier.commonConfig();
            final String outputDirectory = config.outputDirectory();
            final String tumorSample = config.tumorSample();

            // JOBA: Load BAFs
            final BAFSupplier bafSupplier = new BAFSupplier(config, configSupplier.bafConfig());
            final Multimap<String, TumorBAF> bafs = bafSupplier.get();
            final Gender gender = Gender.fromBAFCount(bafs);
            LOGGER.info("Sample gender is {}", gender.toString().toLowerCase());

            // JOBA: Load structural and somatic variants
            final List<SomaticVariant> somaticVariants = somaticVariants(configSupplier);
            final List<StructuralVariant> structuralVariants = structuralVariants(configSupplier);

            final RatioSupplier ratioSupplier;
            final List<GenomeRegion> regions;
            if (cmd.hasOption(COBALT)) {
                final Multimap<String, GCContent> gcContent = FreecGCContentFactory.loadGCContent(cmd.getOptionValue(GC_PROFILE));
                ratioSupplier = new ReadCountRatioSupplier(config, gcContent);
                final Map<String, ChromosomeLength> lengths = new ChromosomeLengthSupplier(config, ratioSupplier.tumorRatios()).get();
                regions = new PCFSegmentSupplier(executorService, config, lengths).get();
            } else {
                final FreecRatioSupplier freecRatioSupplier = new FreecRatioSupplier(config);
                ratioSupplier = freecRatioSupplier;
                regions = new FreecSegmentSupplier(freecRatioSupplier).get();
            }

            LOGGER.info("Merging structural variants into freec segmentation");
            final List<PurpleSegment> segments = PurpleSegmentFactory.createSegments(regions, structuralVariants);

            LOGGER.info("Mapping all observations to the segmented regions");
            final ObservedRegionFactory observedRegionFactory = new ObservedRegionFactory(gender);
            final List<ObservedRegion> observedRegions =
                    observedRegionFactory.combine(segments, bafs, ratioSupplier.tumorRatios(), ratioSupplier.referenceRatios());

            final double cnvRatioWeight = defaultValue(cmd, CNV_RATIO_WEIGHT_FACTOR, CNV_RATIO_WEIGHT_FACTOR_DEFAULT);
            final boolean ploidyPenaltyExperiment = cmd.hasOption(PLOIDY_PENALTY_EXPERIMENT);
            final double observedBafExponent = defaultValue(cmd, OBSERVED_BAF_EXPONENT, OBSERVED_BAF_EXPONENT_DEFAULT);
            final FittedRegionFactory fittedRegionFactory =
                    new FittedRegionFactory(gender, MAX_PLOIDY, cnvRatioWeight, ploidyPenaltyExperiment, observedBafExponent);

            LOGGER.info("Fitting purity");
            final double minPurity = defaultValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
            final double maxPurity = defaultValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
            final double minNormFactor = defaultValue(cmd, MIN_NORM_FACTOR, MIN_NORM_FACTOR_DEFAULT);
            final double maxNormFactor = defaultValue(cmd, MAX_NORM_FACTOR, MAX_NORM_FACTOR_DEFAULT);
            final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(executorService,
                    MAX_PLOIDY,
                    minPurity,
                    maxPurity,
                    PURITY_INCREMENTS,
                    minNormFactor,
                    maxNormFactor,
                    NORM_FACTOR_INCREMENTS,
                    fittedRegionFactory,
                    observedRegions);

            final BestFitFactory bestFitFactory = new BestFitFactory(fittedPurityFactory.bestFitPerPurity(), somaticVariants);
            final FittedPurity bestFit = bestFitFactory.bestFit();
            final List<FittedRegion> fittedRegions = fittedRegionFactory.fitRegion(bestFit.purity(), bestFit.normFactor(), observedRegions);

            final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, bestFit.purity(), bestFit.normFactor());
            final PurpleCopyNumberFactory purpleCopyNumberFactory = new PurpleCopyNumberFactory(purityAdjuster, fittedRegions);
            final List<PurpleCopyNumber> highConfidence = purpleCopyNumberFactory.highConfidenceRegions();
            final List<PurpleCopyNumber> smoothRegions = purpleCopyNumberFactory.smoothedRegions();
            final List<GeneCopyNumber> geneCopyNumbers = geneCopyNumbers(smoothRegions);
            final List<FittedRegion> enrichedFittedRegions = updateRegionsWithCopyNumbers(fittedRegions, highConfidence, smoothRegions);

            final PurityContext purityContext = ImmutablePurityContext.builder()
                    .bestFit(bestFitFactory.bestFit())
                    .bestPerPurity(fittedPurityFactory.bestFitPerPurity())
                    .status(bestFitFactory.status())
                    .gender(gender)
                    .score(bestFitFactory.score())
                    .polyClonalProportion(polyclonalProproption(smoothRegions))
                    .build();

            if (cmd.hasOption(DB_ENABLED)) {
                LOGGER.info("Persisting to database");
                final DatabaseAccess dbAccess = databaseAccess(cmd);
                dbAccess.writePurity(tumorSample, purityContext);
                dbAccess.writeCopynumbers(tumorSample, smoothRegions);
                dbAccess.writeCopynumberRegions(tumorSample, enrichedFittedRegions);
                dbAccess.writeGeneCopynumberRegions(tumorSample, geneCopyNumbers);
                dbAccess.writeStructuralVariants(tumorSample, structuralVariants);
            }

            LOGGER.info("Writing purple data to: {}", outputDirectory);
            FittedPurityFile.write(outputDirectory, tumorSample, purityContext);
            PurpleCopyNumberFile.write(outputDirectory, tumorSample, smoothRegions);
            FittedRegionWriter.writeCopyNumber(outputDirectory, tumorSample, enrichedFittedRegions);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilename(outputDirectory, tumorSample), geneCopyNumbers);

            if (cmd.hasOption(PLOT)) {
                final String baseCircosDirectory = config.outputDirectory() + File.separator + "circos";
                LOGGER.info("Writing circos data to: {}", baseCircosDirectory);
                final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory =
                        new EnrichedSomaticVariantFactory(purityContext.bestFit(), smoothRegions);
                final List<EnrichedSomaticVariant> enrichedSomatics = enrichedSomaticVariantFactory.enrich(somaticVariants);

                new GenerateCircosDataHelper(tumorSample, baseCircosDirectory).write(gender,
                        smoothRegions,
                        enrichedSomatics,
                        structuralVariants);

                final String basePlotDirectory = config.outputDirectory() + File.separator + "plot";
                LOGGER.info("Writing plots to: {}", basePlotDirectory);
                new ChartWriter(tumorSample, basePlotDirectory).write(purityContext.bestFit(), purityContext.score(), smoothRegions, enrichedSomatics);
            }

        } finally {
            executorService.shutdown();
        }
        LOGGER.info("Complete");
    }

    @NotNull
    private List<StructuralVariant> structuralVariants(@NotNull final ConfigSupplier configSupplier) throws IOException {
        final StructuralVariantConfig config = configSupplier.structuralVariantConfig();
        if (config.file().isPresent()) {
            final String filePath = config.file().get().toString();
            LOGGER.info("Loading structural variants from {}", filePath);
            return StructuralVariantFileLoader.fromFile(filePath);
        } else {
            LOGGER.info("Structural variants support disabled.");
            return Collections.emptyList();
        }
    }

    @NotNull
    private List<SomaticVariant> somaticVariants(@NotNull final ConfigSupplier configSupplier) throws IOException, HartwigException {
        final SomaticConfig config = configSupplier.somaticConfig();
        if (config.file().isPresent()) {
            String filename = config.file().get().toString();
            LOGGER.info("Loading somatic variants from {}", filename);
            return VCFFileLoader.loadSomaticVCF(filename)
                    .variants()
                    .stream()
                    .filter(x -> x.type() == VariantType.SNP)
                    .filter(VariantFilter::isPass)
                    .collect(Collectors.toList());
        } else {
            LOGGER.info("Somatic variants support disabled.");
            return Collections.emptyList();
        }
    }

    private List<GeneCopyNumber> geneCopyNumbers(List<PurpleCopyNumber> copyNumbers) throws IOException, EmptyFileException {
        final InputStream inputStream = getClass().getResourceAsStream("/bed/hmf_gene_panel.tsv");
        final List<HmfGenomeRegion> hmgGenomeRegions = HmfSlicerFileLoader.fromInputStream(inputStream);
        return GeneCopyNumberFactory.geneCopyNumbers(hmgGenomeRegions, copyNumbers);
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
    private static Options createOptions() {
        final Options options = new Options();
        ConfigSupplier.addOptions(options);

        options.addOption(OBSERVED_BAF_EXPONENT, true, "Observed baf exponent. Default 1");
        options.addOption(PLOIDY_PENALTY_EXPERIMENT, false, "Use experimental ploidy penality.");
        options.addOption(CNV_RATIO_WEIGHT_FACTOR, true, "CNV ratio deviation scaling.");

        options.addOption(MIN_PURITY, true, "Minimum purity (default 0.05)");
        options.addOption(MAX_PURITY, true, "Maximum purity (default 1.0)");

        options.addOption(MIN_NORM_FACTOR, true, "Minimum norm factor (default 0.33)");
        options.addOption(MAX_NORM_FACTOR, true, "Maximum norm factor (default 2.0)");

        options.addOption(DB_ENABLED, false, "Persist data to DB.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(THREADS, true, "Number of threads (default 2)");
        options.addOption(COBALT, false, "Use cobalt segmentation.");
        options.addOption(GC_PROFILE, true, "Location of GC Profile.");
        options.addOption(PLOT, false, "Generate circos data.");

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
