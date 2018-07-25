package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFactory.polyclonalProportion;
import static com.hartwig.hmftools.patientdb.LoadPurpleData.persistToDatabase;
import static com.hartwig.hmftools.purple.PurpleRegionZipper.updateRegionsWithCopyNumbers;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.gc.GCProfileFactory;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFactory;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.BestFitFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFactory;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactoryV1;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactoryV2;
import com.hartwig.hmftools.common.purple.region.FittedRegionFile;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionFactory;
import com.hartwig.hmftools.common.purple.segment.Cluster;
import com.hartwig.hmftools.common.purple.segment.ClusterFactory;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.purple.segment.PurpleSegmentFactory;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.version.VersionInfo;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.config.CircosConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.DBConfig;
import com.hartwig.hmftools.purple.config.FittingConfig;
import com.hartwig.hmftools.purple.config.SmoothingConfig;
import com.hartwig.hmftools.purple.config.SomaticConfig;
import com.hartwig.hmftools.purple.config.StructuralVariantConfig;
import com.hartwig.hmftools.purple.plot.ChartWriter;
import com.hartwig.hmftools.purple.ratio.ChromosomeLengthSupplier;
import com.hartwig.hmftools.purple.segment.PCFPositionsSupplier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class PurityPloidyEstimateApplication {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    private static final int THREADS_DEFAULT = 2;

    private static final String THREADS = "threads";
    private static final String EXPERIMENTAL = "experimental";
    private static final String PLOIDY_PENALTY_FACTOR = "ploidy_penalty_factor";
    private static final String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    private static final String VERSION = "version";

    private static final String CNV_RATIO_WEIGHT_FACTOR = "cnv_ratio_weight_factor";
    private static final double CNV_RATIO_WEIGHT_FACTOR_DEFAULT = 0.2;

    private static final String OBSERVED_BAF_EXPONENT = "observed_baf_exponent";
    private static final double OBSERVED_BAF_EXPONENT_DEFAULT = 1;

    public static void main(final String... args)
            throws ParseException, IOException, SQLException, ExecutionException, InterruptedException {
        new PurityPloidyEstimateApplication(args);
    }

    private PurityPloidyEstimateApplication(final String... args)
            throws ParseException, IOException, SQLException, ExecutionException, InterruptedException {
        final VersionInfo version = new VersionInfo("purple.version");
        LOGGER.info("PURPLE version: {}", version.version());

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (cmd.hasOption(VERSION)) {
            System.exit(0);
        }

        final int threads = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : THREADS_DEFAULT;
        final ExecutorService executorService = Executors.newFixedThreadPool(threads);
        try {
            // JOBA: Get common config
            final ConfigSupplier configSupplier = new ConfigSupplier(cmd, options);
            final CommonConfig config = configSupplier.commonConfig();
            final String outputDirectory = config.outputDirectory();
            final String tumorSample = config.tumorSample();

            // JOBA: Read Gene Panel
            final List<HmfGenomeRegion> genePanel = HmfGenePanelSupplier.allGeneList();

            // JOBA: Load BAFs from AMBER
            final String amberFile = configSupplier.bafConfig().bafFile().toString();
            LOGGER.info("Reading amber bafs from {}", amberFile);
            final Multimap<String, AmberBAF> bafs = AmberBAFFile.read(amberFile);
            int averageTumorDepth =
                    (int) Math.round(bafs.values().stream().mapToInt(AmberBAF::tumorDepth).filter(x -> x > 0).average().orElse(100));
            LOGGER.info("Average amber tumor depth is {} reads", averageTumorDepth);

            // JOBA: Load Ratios from COBALT
            final String ratioFilename = CobaltRatioFile.generateFilename(config.cobaltDirectory(), config.tumorSample());
            LOGGER.info("Reading cobalt ratios from {}", ratioFilename);
            final ListMultimap<String, CobaltRatio> ratios = CobaltRatioFile.read(ratioFilename);

            LOGGER.info("Reading GC Profiles from {}", config.gcProfile());
            final Multimap<String, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(config.gcProfile());

            // JOBA: Gender
            final Gender amberGender = Gender.fromAmber(bafs);
            final Gender cobaltGender = Gender.fromCobalt(ratios);
            if (cobaltGender.equals(amberGender)) {
                LOGGER.info("Sample gender is {}", cobaltGender.toString().toLowerCase());
            } else {
                LOGGER.warn("COBALT gender {} does not match AMBER gender {}", cobaltGender, amberGender);
            }

            // JOBA: Load structural and somatic variants
            final List<SomaticVariant> somaticVariants = somaticVariants(configSupplier);
            final List<StructuralVariant> structuralVariants = structuralVariants(configSupplier);

            // JOBA: Ratio Segmentation
            final Map<String, ChromosomeLength> lengths = new ChromosomeLengthSupplier(config, ratios).get();
            final Multimap<String, PCFPosition> pcfPositions = PCFPositionsSupplier.createPositions(config);
            final Multimap<String, Cluster> clusterMap =
                    new ClusterFactory(config.windowSize()).cluster(structuralVariants, pcfPositions, ratios);
            final List<PurpleSegment> segments = PurpleSegmentFactory.segment(clusterMap, lengths);

            LOGGER.info("Mapping all observations to the segmented regions");
            final ObservedRegionFactory observedRegionFactory = new ObservedRegionFactory(config.windowSize(), cobaltGender);
            final List<ObservedRegion> observedRegions = observedRegionFactory.combine(segments, bafs, ratios, gcProfiles);

            LOGGER.info("Fitting purity");
            final FittingConfig fittingConfig = configSupplier.fittingConfig();
            final double cnvRatioWeight = defaultValue(cmd, CNV_RATIO_WEIGHT_FACTOR, CNV_RATIO_WEIGHT_FACTOR_DEFAULT);
            final double observedBafExponent = defaultValue(cmd, OBSERVED_BAF_EXPONENT, OBSERVED_BAF_EXPONENT_DEFAULT);

            final double ploidyPenaltyFactor = defaultValue(cmd, PLOIDY_PENALTY_FACTOR, 0.25);
            final double ploidyPenaltyStandardDevation = defaultValue(cmd, PLOIDY_PENALTY_STANDARD_DEVIATION, 0.03);

            final FittedRegionFactory fittedRegionFactory = cmd.hasOption(EXPERIMENTAL)
                    ? new FittedRegionFactoryV2(cobaltGender, averageTumorDepth, ploidyPenaltyFactor, ploidyPenaltyStandardDevation)
                    : new FittedRegionFactoryV1(cobaltGender,
                            fittingConfig.maxPloidy(),
                            cnvRatioWeight,
                            averageTumorDepth,
                            observedBafExponent);

            final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(executorService,
                    fittingConfig.maxPloidy(),
                    fittingConfig.minPurity(),
                    fittingConfig.maxPurity(),
                    fittingConfig.purityIncrement(),
                    fittingConfig.minNormFactor(),
                    fittingConfig.maxNormFactor(),
                    fittingConfig.normFactorIncrement(),
                    fittedRegionFactory,
                    observedRegions);

            final List<FittedPurity> bestFitPerPurity = fittedPurityFactory.bestFitPerPurity();

            final SomaticConfig somaticConfig = configSupplier.somaticConfig();
            final List<SomaticVariant> snps =
                    somaticVariants.stream().filter(x -> x.type() == VariantType.SNP).collect(Collectors.toList());
            final BestFitFactory bestFitFactory =
                    new BestFitFactory(somaticConfig.minTotalVariants(), somaticConfig.minPeakVariants(), bestFitPerPurity, snps);
            final FittedPurity bestFit = bestFitFactory.bestFit();

            final List<FittedRegion> fittedRegions = fittedRegionFactory.fitRegion(bestFit.purity(), bestFit.normFactor(), observedRegions);

            final PurityAdjuster purityAdjuster = new PurityAdjuster(cobaltGender, bestFit.purity(), bestFit.normFactor());

            final SmoothingConfig smoothingConfig = configSupplier.smoothingConfig();
            final PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(smoothingConfig.minDiploidTumorRatioCount(),
                    smoothingConfig.minDiploidTumorRatioCountAtCentromere(),
                    cobaltGender,
                    purityAdjuster,
                    fittedRegions,
                    structuralVariants);
            final List<PurpleCopyNumber> copyNumbers = copyNumberFactory.copyNumbers();
            final List<PurpleCopyNumber> germlineDeletions = copyNumberFactory.germlineDeletions();

            final List<FittedRegion> enrichedFittedRegions = updateRegionsWithCopyNumbers(fittedRegions, copyNumbers);

            final PurityContext purityContext = ImmutablePurityContext.builder()
                    .version(version.version())
                    .bestFit(bestFitFactory.bestFit())
                    .status(bestFitFactory.status())
                    .gender(cobaltGender)
                    .score(bestFitFactory.score())
                    .polyClonalProportion(polyclonalProportion(copyNumbers))
                    .build();

            final List<PurityAdjustedSomaticVariant> enrichedSomatics =
                    new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, enrichedFittedRegions).create(somaticVariants);

            final List<GeneCopyNumber> geneCopyNumbers =
                    GeneCopyNumberFactory.geneCopyNumbers(genePanel, copyNumbers, germlineDeletions, enrichedSomatics);

            LOGGER.info("Generating QC Stats");
            final PurpleQC qcChecks =
                    PurpleQCFactory.create(bestFitFactory.bestFit(), copyNumbers, amberGender, cobaltGender, geneCopyNumbers);

            final DBConfig dbConfig = configSupplier.dbConfig();
            if (dbConfig.enabled()) {
                final DatabaseAccess dbAccess = databaseAccess(dbConfig);
                persistToDatabase(dbAccess,
                        tumorSample,
                        bestFitPerPurity,
                        copyNumbers,
                        germlineDeletions,
                        enrichedFittedRegions,
                        purityContext,
                        qcChecks,
                        geneCopyNumbers);
            }

            LOGGER.info("Writing purple data to: {}", outputDirectory);
            version.write(outputDirectory);
            PurpleQCFile.write(PurpleQCFile.generateFilename(outputDirectory, tumorSample), qcChecks);
            FittedPurityFile.write(outputDirectory, tumorSample, purityContext);
            FittedPurityRangeFile.write(outputDirectory, tumorSample, bestFitPerPurity);
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilename(outputDirectory, tumorSample), copyNumbers);
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateGermlineFilename(outputDirectory, tumorSample), germlineDeletions);
            FittedRegionFile.write(FittedRegionFile.generateFilename(outputDirectory, tumorSample), enrichedFittedRegions);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilename(outputDirectory, tumorSample), geneCopyNumbers);

            final CircosConfig circosConfig = configSupplier.circosConfig();
            LOGGER.info("Writing plots to: {}", circosConfig.plotDirectory());
            new ChartWriter(tumorSample, circosConfig.plotDirectory()).write(purityContext.bestFit(),
                    purityContext.score(),
                    copyNumbers,
                    enrichedSomatics);

            LOGGER.info("Writing circos data to: {}", circosConfig.circosDirectory());
            new GenerateCircosData(configSupplier, executorService).write(cobaltGender,
                    copyNumbers,
                    enrichedSomatics,
                    structuralVariants,
                    fittedRegions,
                    Lists.newArrayList(bafs.values()));
        } finally {
            executorService.shutdown();
        }
        LOGGER.info("Complete");
    }

    @NotNull
    private static List<StructuralVariant> structuralVariants(@NotNull final ConfigSupplier configSupplier) throws IOException {
        final StructuralVariantConfig config = configSupplier.structuralVariantConfig();
        if (config.file().isPresent()) {
            final String filePath = config.file().get().toString();
            LOGGER.info("Loading structural variants from {}", filePath);
            return StructuralVariantFileLoader.fromFile(filePath, true);
        } else {
            LOGGER.info("Structural variants support disabled.");
            return Collections.emptyList();
        }
    }

    @NotNull
    private static List<SomaticVariant> somaticVariants(@NotNull final ConfigSupplier configSupplier) throws IOException {
        final SomaticConfig config = configSupplier.somaticConfig();
        if (config.file().isPresent()) {
            String filename = config.file().get().toString();
            LOGGER.info("Loading somatic variants from {}", filename);

            SomaticVariantFactory factory =
                    SomaticVariantFactory.filteredInstance(new PassingVariantFilter(), new NTFilter(), new SGTFilter());

            return factory.fromVCFFile(configSupplier.commonConfig().tumorSample(), filename);
        } else {
            LOGGER.info("Somatic variants support disabled.");
            return Collections.emptyList();
        }
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
        options.addOption(CNV_RATIO_WEIGHT_FACTOR, true, "CNV ratio deviation scaling.");

        options.addOption(THREADS, true, "Number of threads (default 2)");
        options.addOption(EXPERIMENTAL, false, "Anything goes!");
        options.addOption(VERSION, false, "Exit after displaying version info.");
        options.addOption(PLOIDY_PENALTY_FACTOR, true, "Ploidy Penalty Factor");
        options.addOption(PLOIDY_PENALTY_STANDARD_DEVIATION, true, "Ploidy Penalty Standard Deviation");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final DBConfig dbConfig) throws SQLException {
        return new DatabaseAccess(dbConfig.user(), dbConfig.password(), dbConfig.url());
    }
}
