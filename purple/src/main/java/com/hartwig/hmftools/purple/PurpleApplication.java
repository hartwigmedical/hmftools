package com.hartwig.hmftools.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSvFile;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MAX_DELETED_GENES;
import static com.hartwig.hmftools.common.purple.GeneCopyNumber.listToMap;
import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NORMAL;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.MemoryCalcs.calcMemoryUsage;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.hasDatabaseConfig;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleSummaryData.createPurity;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_MAX_DELETED_GENES;
import static com.hartwig.hmftools.purple.gene.PurpleRegionZipper.updateRegionsWithCopyNumbers;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableBestFit;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.purple.gene.GeneCopyNumberBuilder;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.purple.purity.PurityAdjusterAbnormalChromosome;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.BestFit;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.SegmentFile;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.FittingConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.config.SampleData;
import com.hartwig.hmftools.purple.config.SampleDataFiles;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;
import com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.purple.germline.GermlineDeletionDrivers;
import com.hartwig.hmftools.purple.germline.GermlineDrivers;
import com.hartwig.hmftools.purple.fitting.BestFitFactory;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.purple.fitting.PeakModelFile;
import com.hartwig.hmftools.purple.germline.GermlineVariants;
import com.hartwig.hmftools.purple.plot.Charts;
import com.hartwig.hmftools.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.purple.recovery.RecoverStructuralVariants;
import com.hartwig.hmftools.purple.somatic.SomaticPeakStream;
import com.hartwig.hmftools.purple.somatic.SomaticPurityEnrichment;
import com.hartwig.hmftools.purple.somatic.SomaticStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class PurpleApplication
{
    private final VersionInfo mPurpleVersion;
    private final ExecutorService mExecutorService;
    private final ReferenceData mReferenceData;
    private final CommandLine mCmdLineArgs;
    private final PurpleConfig mConfig;

    private final GermlineVariants mGermlineVariants;
    private final Segmentation mSegmentation;
    private final Charts mCharts;

    private final DatabaseAccess mDbAccess;

    private static final int THREADS_DEFAULT = 2;
    private static final String THREADS = "threads";
    private static final String VERSION = "version";

    private PurpleApplication(final Options options, final String... args) throws ParseException, IOException
    {
        mPurpleVersion = new VersionInfo("purple.version");
        PPL_LOGGER.info("Purple version: {}", mPurpleVersion.version());

        mCmdLineArgs = createCommandLine(options, args);

        if(mCmdLineArgs.hasOption(VERSION))
        {
            System.exit(0);
        }

        setLogLevel(mCmdLineArgs);

        mDbAccess = hasDatabaseConfig(mCmdLineArgs) ? createDatabaseAccess(mCmdLineArgs) : null;

        // load config
        mConfig = new PurpleConfig(mPurpleVersion.version(), mCmdLineArgs);

        if(!mConfig.isValid())
        {
            PPL_LOGGER.error("initialisation error, exiting");
            System.exit(1);
        }

        // and common reference data
        mReferenceData = new ReferenceData(mCmdLineArgs, mConfig);

        if(!mReferenceData.isValid())
        {
            PPL_LOGGER.error("initialisation error, exiting");
            System.exit(1);
        }

        final int threads = mCmdLineArgs.hasOption(THREADS) ? Integer.parseInt(mCmdLineArgs.getOptionValue(THREADS)) : THREADS_DEFAULT;
        mExecutorService = Executors.newFixedThreadPool(threads);

        mGermlineVariants = new GermlineVariants(mConfig, mReferenceData, mPurpleVersion.version());

        if(!mConfig.DriversOnly)
        {
            mSegmentation = new Segmentation(mReferenceData);
            mCharts = new Charts(mConfig, mExecutorService, mReferenceData.RefGenVersion.is38());
        }
        else
        {
            mSegmentation = null;
            mCharts = null;
        }
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        try
        {
            processSample(mConfig.ReferenceId, mConfig.TumorId);
        }
        finally
        {
            mExecutorService.shutdown();
        }

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        PPL_LOGGER.info("Purple complete, runTime({})", format("%.1fs", timeTakenMs/1000.0));
    }

    private SampleData loadSampleData(final String referenceId, final String tumorId, final SampleDataFiles sampleDataFiles) throws Exception
    {
        SampleData sampleData = null;

        final SomaticVariantCache somaticVariantCache = new SomaticVariantCache(mConfig);

        if(!mConfig.DriversOnly)
        {
            // load amber and cobalt sample data
            final AmberData amberData = new AmberData(
                    mConfig.germlineMode() ? referenceId : tumorId, sampleDataFiles.AmberDirectory, mConfig.germlineMode(),
                    mReferenceData.RefGenVersion);

            final CobaltData cobaltData = new CobaltData(
                    referenceId, tumorId, sampleDataFiles.CobaltDirectory, amberData.PatientGender,
                    mConfig.tumorOnlyMode(), mConfig.germlineMode());

            // load structural and somatic variants
            final StructuralVariantCache svCache = createStructuralVariantCache(tumorId, sampleDataFiles);

            sampleData = new SampleData(referenceId, tumorId, amberData, cobaltData, svCache, somaticVariantCache);
        }
        else
        {
            sampleData = new SampleData(referenceId, tumorId, null, null, null, somaticVariantCache);
        }

        if(mConfig.runTumor())
            somaticVariantCache.loadSomatics(sampleDataFiles.SomaticVcfFile, mReferenceData.SomaticHotspots);

        return sampleData;
    }

    private void processSample(final String referenceId, final String tumorSample)
    {
        try
        {
            if(mConfig.DriversOnly)
            {
                findDrivers(tumorSample);
            }
            else
            {
                final SampleDataFiles sampleDataFiles = new SampleDataFiles(mCmdLineArgs, tumorSample);
                final SampleData sampleData = loadSampleData(referenceId, tumorSample, sampleDataFiles);

                if(sampleData == null)
                    System.exit(1);

                PPL_LOGGER.debug("post-data-loading memory({}mb)", calcMemoryUsage());

                performFit(referenceId, tumorSample, sampleDataFiles, sampleData);
            }
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed processing sample({}): {}", tumorSample, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void performFit(
            final String referenceId, final String tumorSample,
            final SampleDataFiles sampleDataFiles, final SampleData sampleData) throws Exception
    {
        final AmberData amberData = sampleData.Amber;
        final CobaltData cobaltData = sampleData.Cobalt;

        final Gender amberGender = amberData.PatientGender;
        final Gender cobaltGender = cobaltData.gender();
        final Gender gender = mConfig.tumorOnlyMode() ? amberGender : cobaltData.gender();

        if(cobaltGender.equals(amberGender))
        {
            PPL_LOGGER.info("sample gender is {}", cobaltGender.toString().toLowerCase());
        }
        else
        {
            PPL_LOGGER.warn("Cobalt gender {} does not match Amber gender {}", cobaltGender, amberGender);
        }

        final CobaltChromosomes cobaltChromosomes = cobaltData.CobaltChromosomes;
        final SomaticVariantCache somaticCache = mConfig.runTumor() ? sampleData.SomaticCache : null;

        PPL_LOGGER.info("applying segmentation");
        final List<ObservedRegion> observedRegions = mSegmentation.createSegments(sampleData.SvCache.variants(), amberData, cobaltData);

        if(observedRegions.isEmpty())
        {
            PPL_LOGGER.warn("no observed regions created");
            System.exit(0);
        }

        PPL_LOGGER.info("purple output directory: {}", mConfig.OutputDir);
        mPurpleVersion.write(mConfig.OutputDir);

        final List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        final List<ObservedRegion> fittedRegions = Lists.newArrayList();

        BestFit bestFit = null;
        FittedPurity fittedPurity = null;
        PurityAdjuster purityAdjuster = null;
        Set<String> reportedGenes = Sets.newHashSet();
        SomaticStream somaticStream = null;

        final FittedRegionFactory fittedRegionFactory =
                createFittedRegionFactory(amberData.AverageTumorDepth, cobaltChromosomes, mConfig.Fitting);

        if(mConfig.runTumor())
        {
            PPL_LOGGER.info("fitting purity");

            bestFit = fitPurity(sampleData, observedRegions, fittedRegionFactory, sampleData.SvCache.variants());

            fittedPurity = bestFit.fit();

            purityAdjuster =
                    new PurityAdjusterAbnormalChromosome(fittedPurity.purity(), fittedPurity.normFactor(), cobaltChromosomes.chromosomes());

            final PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(
                    mConfig.Fitting.MinDiploidTumorRatioCount,
                    mConfig.Fitting.MinDiploidTumorRatioCountAtCentromere,
                    amberData.AverageTumorDepth,
                    fittedPurity.ploidy(),
                    purityAdjuster,
                    cobaltData.CobaltChromosomes);

            PPL_LOGGER.info("calculating copy number");
            fittedRegions.addAll(fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), observedRegions));

            copyNumberFactory.invoke(fittedRegions, sampleData.SvCache.variants());

            final int recoveredSVCount = recoverStructuralVariants(
                    sampleData, sampleDataFiles, purityAdjuster, copyNumberFactory.copyNumbers());

            if(recoveredSVCount > 0)
            {
                PPL_LOGGER.info("reapplying segmentation with {} recovered structural variants", recoveredSVCount);
                final List<ObservedRegion> recoveredObservedRegions =
                        mSegmentation.createSegments(sampleData.SvCache.variants(), amberData, cobaltData);

                PPL_LOGGER.info("recalculating copy number");
                fittedRegions.clear();
                fittedRegions.addAll(fittedRegionFactory.fitRegion(
                        fittedPurity.purity(), fittedPurity.normFactor(), recoveredObservedRegions));

                copyNumberFactory.invoke(fittedRegions, sampleData.SvCache.variants());
            }

            copyNumbers.addAll(copyNumberFactory.copyNumbers());
            sampleData.SvCache.inferMissingVariant(copyNumbers);

            final List<ObservedRegion> enrichedObservedRegions = updateRegionsWithCopyNumbers(fittedRegions, copyNumbers);

            geneCopyNumbers.addAll(GeneCopyNumberBuilder.createGeneCopyNumbers(
                    mReferenceData.RefGenVersion, mReferenceData.GeneTransCache, copyNumbers));

            PPL_LOGGER.debug("post-fit memory({}mb)", calcMemoryUsage());

            final List<PeakModel> somaticPeaks = Lists.newArrayList();

            PPL_LOGGER.info("modelling somatic peaks");
            final SomaticPeakStream somaticPeakStream = new SomaticPeakStream(mConfig);

            final SomaticPurityEnrichment somaticPurityEnrichment = new SomaticPurityEnrichment(
                    mConfig.Version, mConfig.TumorId, purityAdjuster, copyNumbers, fittedRegions);

            sampleData.SomaticCache.purityEnrich(somaticPurityEnrichment);

            List<PeakModel> peakModelValues = somaticPeakStream.somaticPeakModel(somaticCache);
            somaticPeaks.addAll(peakModelValues);

            // at the moment the enriching of somatic variants is also contributing to the purity context, so it cannot be done afterwards
            // if the read and write process were split then so could the fitting and enriching steps
            PPL_LOGGER.info("enriching somatic variants");

            somaticStream = new SomaticStream(mConfig, mReferenceData, somaticCache, somaticPeaks);

            somaticStream.processAndWrite(purityAdjuster, copyNumbers, enrichedObservedRegions);

            PPL_LOGGER.debug("post-enrichment memory({}mb)", calcMemoryUsage());

            sampleData.SvCache.write(purityAdjuster, copyNumbers, mConfig.tumorOnlyMode());

            reportedGenes.addAll(somaticStream.reportedGenes());

            FittedPurityRangeFile.write(mConfig.OutputDir, tumorSample, bestFit.allFits());
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample), copyNumbers);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample), geneCopyNumbers);
            PeakModelFile.write(PeakModelFile.generateFilename(mConfig.OutputDir, tumorSample), somaticPeaks);
        }
        else
        {
            fittedPurity = ImmutableFittedPurity.builder()
                    .purity(1)
                    .ploidy(2)
                    .normFactor(1)
                    .diploidProportion(1)
                    .somaticPenalty(0)
                    .score(0)
                    .build();

            FittedPurityScore score = ImmutableFittedPurityScore.builder()
                    .minPloidy(fittedPurity.ploidy())
                    .maxPloidy(fittedPurity.ploidy())
                    .minPurity(fittedPurity.purity())
                    .maxPurity(fittedPurity.purity())
                    .minDiploidProportion(1)
                    .maxDiploidProportion(1)
                    .build();

            bestFit = ImmutableBestFit.builder()
                    .fit(fittedPurity)
                    .method(NORMAL)
                    .score(score)
                    .allFits(Lists.newArrayList(fittedPurity))
                    .build();

            fittedRegions.addAll(fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), observedRegions));
        }

        PPL_LOGGER.info("generating QC Stats");
        final PurpleQC qcChecks = PurpleSummaryData.createQC(
                amberData.Contamination, bestFit, amberGender, cobaltGender, copyNumbers, geneCopyNumbers,
                cobaltChromosomes.germlineAberrations(), amberData.AverageTumorDepth,
                mConfig.TargetRegionsMode ? TARGET_REGIONS_MAX_DELETED_GENES : MAX_DELETED_GENES);

        final PurityContext purityContext = createPurity(
                mPurpleVersion.version(), bestFit, gender, mConfig, qcChecks, copyNumbers, somaticStream, sampleData.SvCache);

        PurityContextFile.write(mConfig.OutputDir, tumorSample, purityContext);
        SegmentFile.write(SegmentFile.generateFilename(mConfig.OutputDir, tumorSample), fittedRegions);

        List<GermlineDeletion> germlineDeletions = Lists.newArrayList();
        if(mConfig.runGermline())
        {
            mGermlineVariants.processAndWrite(
                    referenceId, tumorSample, sampleDataFiles.GermlineVcfFile, purityAdjuster, copyNumbers, reportedGenes);

            GermlineDeletionDrivers germlineDeletionDrivers = new GermlineDeletionDrivers(
                    mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

            germlineDeletionDrivers.findDeletions(copyNumbers, fittedRegions);
            germlineDeletions.addAll(germlineDeletionDrivers.getDeletions());
            GermlineDeletion.write(GermlineDeletion.generateFilename(mConfig.OutputDir, tumorSample), germlineDeletions);
        }

        if(mDbAccess != null)
        {
            PPL_LOGGER.info("writing purple data to database: {}", mCmdLineArgs.getOptionValue(DB_URL));

            mDbAccess.writePurity(tumorSample, purityContext, qcChecks);

            if(mConfig.runTumor())
            {
                mDbAccess.writeBestFitPerPurity(tumorSample, bestFit.bestFitPerPurity());
                mDbAccess.writeCopynumbers(tumorSample, copyNumbers);
                mDbAccess.writeGeneCopyNumbers(tumorSample, geneCopyNumbers);
            }

            if(mConfig.runGermline())
                mDbAccess.writeGermlineDeletions(tumorSample, germlineDeletions);
        }

        if(!mConfig.germlineMode() && (mConfig.Charting.Enabled || mConfig.Charting.CircosBinary.isPresent()))
        {
            PPL_LOGGER.info("generating charts");

            mCharts.write(
                    referenceId, tumorSample, !sampleDataFiles.SomaticVcfFile.isEmpty(),
                    gender, copyNumbers, somaticStream.downsampledVariants(), sampleData.SvCache.variants(),
                    fittedRegions, Lists.newArrayList(amberData.ChromosomeBafs.values()));
        }

        if(mConfig.RunDrivers)
        {
            findDrivers(tumorSample, purityContext, copyNumbers, geneCopyNumbers, fittedRegions, somaticStream);
        }
    }

    private void findDrivers(final String tumorSample) throws Exception
    {
        final String purpleDataPath = mConfig.OutputDir;
        SomaticStream somaticStream = null;

        if(mConfig.runTumor())
        {
            String somaticVcf = PurpleCommon.purpleSomaticVcfFile(purpleDataPath, tumorSample);

            SomaticVariantCache somaticVariantCache = new SomaticVariantCache(mConfig);
            ListMultimap<Chromosome, VariantHotspot> emptyHotspots = ArrayListMultimap.create(); // already annotated in VCF
            somaticVariantCache.loadSomatics(somaticVcf, emptyHotspots);

            // the counts passed in here are only for down-sampling for charting, which is not relevant for drivers
            somaticStream = new SomaticStream(mConfig, mReferenceData, somaticVariantCache, null);
        }

        if(mConfig.runGermline())
        {
            final String germlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDataPath, tumorSample);

            if(Files.exists(Paths.get(germlineVcf)))
            {
                mGermlineVariants.loadReportableVariants(germlineVcf);
                PPL_LOGGER.info("loaded {} reportable germline variants", mGermlineVariants.reportableVariants().size());
            }
        }

        final PurityContext purityContext = PurityContextFile.read(purpleDataPath, tumorSample);

        final List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                PurpleCopyNumberFile.generateFilenameForReading(purpleDataPath, tumorSample));

        final List<ObservedRegion> fittedRegions = SegmentFile.read(SegmentFile.generateFilename(purpleDataPath, tumorSample)).stream()
                .filter(x -> x.germlineStatus() == HET_DELETION || x.germlineStatus() == HOM_DELETION)
                .collect(Collectors.toList());

        final List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(
                GeneCopyNumberFile.generateFilenameForReading(purpleDataPath, tumorSample));

        findDrivers(tumorSample, purityContext, copyNumbers, geneCopyNumbers, fittedRegions, somaticStream);
    }

    private void findDrivers(
            final String tumorSample, final PurityContext purityContext, final List<PurpleCopyNumber> copyNumbers,
            final List<GeneCopyNumber> geneCopyNumbers, final List<ObservedRegion> fittedRegions,
            @Nullable final SomaticStream somaticStream) throws IOException
    {
        final List<DriverCatalog> somaticDriverCatalog = Lists.newArrayList();
        final List<DriverCatalog> germlineDriverCatalog = Lists.newArrayList();

        PPL_LOGGER.info("generating drivers");

        final Map<String,List<GeneCopyNumber>> geneCopyNumberMap = listToMap(geneCopyNumbers);

        if(mConfig.runTumor())
        {
            somaticDriverCatalog.addAll(somaticStream.buildDrivers(geneCopyNumberMap));

            final CNADrivers cnaDrivers = new CNADrivers(purityContext.qc().status(), mReferenceData.DriverGenes);

            somaticDriverCatalog.addAll(cnaDrivers.deletions(geneCopyNumbers));

            // partial AMPs are only allowed for WGS
            somaticDriverCatalog.addAll(cnaDrivers.amplifications(
                    purityContext.bestFit().ploidy(), geneCopyNumbers, !mConfig.TargetRegionsMode));

            DriverCatalogFile.write(DriverCatalogFile.generateSomaticFilename(mConfig.OutputDir, tumorSample), somaticDriverCatalog);
        }

        if(mConfig.runGermline())
        {
            final GermlineDrivers germlineDrivers = new GermlineDrivers(mReferenceData.DriverGenes.driverGenes());
            germlineDriverCatalog.addAll(germlineDrivers.findDrivers(mGermlineVariants.reportableVariants(), geneCopyNumberMap));

            GermlineDeletionDrivers germlineDeletionDrivers = new GermlineDeletionDrivers(
                    mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

            germlineDeletionDrivers.findDeletions(copyNumbers, fittedRegions);

            germlineDriverCatalog.addAll(germlineDeletionDrivers.getDrivers());

            DriverCatalogFile.write(DriverCatalogFile.generateGermlineFilename(mConfig.OutputDir, tumorSample), germlineDriverCatalog);
        }

        if(mDbAccess != null)
            mDbAccess.writePurpleDriverCatalog(tumorSample, somaticDriverCatalog, germlineDriverCatalog);
    }

    private int recoverStructuralVariants(final SampleData sampleData, final SampleDataFiles sampleDataFiles,
            final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        if(sampleDataFiles.RecoveredSvVcfFile.isEmpty())
            return 0;

        PPL_LOGGER.info("loading recovery candidates from {}", sampleDataFiles.RecoveredSvVcfFile);

        final RecoverStructuralVariants recovery = new RecoverStructuralVariants(
                mConfig.Fitting, purityAdjuster, sampleDataFiles.RecoveredSvVcfFile, copyNumbers);

        final Collection<VariantContext> recoveredVariants = recovery.recoverVariants(sampleData.SvCache.variants());

        if(!recoveredVariants.isEmpty())
        {
            recoveredVariants.forEach(x -> sampleData.SvCache.addVariant(x));
        }

        return recoveredVariants.size();
    }

    private BestFit fitPurity(final SampleData sampleData, final List<ObservedRegion> observedRegions,
            final FittedRegionFactory fittedRegionFactory, final List<StructuralVariant> structuralVariants)
            throws ExecutionException, InterruptedException
    {
        final FittingConfig fittingConfig = mConfig.Fitting;
        final SomaticFitConfig somaticFitConfig = mConfig.SomaticFitting;

        final List<FittedPurity> fitCandidates = Lists.newArrayList();

        final CobaltChromosomes cobaltChromosomes = sampleData.Cobalt.CobaltChromosomes;

        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(
                mExecutorService,
                cobaltChromosomes,
                fittingConfig.MinPurity,
                fittingConfig.MaxPurity,
                fittingConfig.PurityIncrement,
                fittingConfig.MinPloidy,
                fittingConfig.MaxPloidy,
                somaticFitConfig.SomaticPenaltyWeight,
                mConfig.tumorOnlyMode(),
                fittedRegionFactory,
                observedRegions,
                sampleData.SomaticCache.fittingVariants());

        fitCandidates.addAll(fittedPurityFactory.all());

        final BestFitFactory bestFitFactory = new BestFitFactory(
                mConfig,
                sampleData.Amber.minSomaticTotalReadCount(),
                sampleData.Amber.maxSomaticTotalReadCount(),
                fitCandidates,
                sampleData.SomaticCache.fittingVariants(),
                structuralVariants,
                observedRegions);

        return bestFitFactory.bestFit();
    }

    private FittedRegionFactory createFittedRegionFactory(final int averageTumorDepth, final CobaltChromosomes cobaltChromosomes,
            final FittingConfig fitScoreConfig)
    {
        return new FittedRegionFactory(cobaltChromosomes,
                averageTumorDepth,
                fitScoreConfig.PloidyPenaltyFactor,
                fitScoreConfig.PloidyPenaltyStandardDeviation,
                fitScoreConfig.PloidyPenaltyMinStandardDeviationPerPloidy,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneMultiplier,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneAdditional,
                fitScoreConfig.PloidyPenaltyBaselineDeviation);
    }

    @NotNull
    private StructuralVariantCache createStructuralVariantCache(final String tumorSample, final SampleDataFiles sampleDataFiles)
    {
        if(sampleDataFiles.SvVcfFile.isEmpty())
        {
            return new StructuralVariantCache();
        }

        PPL_LOGGER.info("loading structural variants from {}", sampleDataFiles.SvVcfFile);

        final String outputVcf = purpleSvFile(mConfig.OutputDir, tumorSample);

        return new StructuralVariantCache(mPurpleVersion.version(), sampleDataFiles.SvVcfFile, outputVcf, mReferenceData);
    }

    public static void main(final String... args) throws IOException
    {
        final Options options = createOptions();

        try
        {
            PurpleApplication purpleApplication = new PurpleApplication(options, args);
            purpleApplication.run();
        }
        catch(ParseException e)
        {
            PPL_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PurpleApplication", options);
            System.exit(1);
        }
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();
        PurpleConfig.addOptions(options);

        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(THREADS, true, "Number of threads (default 2)");
        options.addOption(VERSION, false, "Exit after displaying version info.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
