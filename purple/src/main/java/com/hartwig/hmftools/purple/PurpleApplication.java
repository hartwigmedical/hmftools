package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_GERMLINE_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.gene.GeneCopyNumber.listToMap;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.hasDatabaseConfig;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.gene.PurpleRegionZipper.updateRegionsWithCopyNumbers;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.wholeGenomeDuplication;
import static com.hartwig.hmftools.purple.purity.FittedPurityScoreFactory.polyclonalProportion;

import java.io.File;
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

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterAbnormalChromosome;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.region.SegmentFile;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.sv.StructuralVariant;
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
import com.hartwig.hmftools.purple.drivers.GermlineDeletionDrivers;
import com.hartwig.hmftools.purple.drivers.GermlineDrivers;
import com.hartwig.hmftools.purple.fitting.BestFitFactory;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.purple.gene.GeneCopyNumberFactory;
import com.hartwig.hmftools.purple.germline.GermlineVariants;
import com.hartwig.hmftools.purple.plot.Charts;
import com.hartwig.hmftools.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.purple.recovery.RecoverStructuralVariants;
import com.hartwig.hmftools.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.purple.region.FittedRegionFactoryV2;
import com.hartwig.hmftools.purple.somatic.SomaticPeakStream;
import com.hartwig.hmftools.purple.somatic.SomaticStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

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
        PPL_LOGGER.info("PURPLE version: {}", mPurpleVersion.version());

        mCmdLineArgs = createCommandLine(options, args);

        if(mCmdLineArgs.hasOption(VERSION))
        {
            System.exit(0);
        }

        setLogLevel(mCmdLineArgs);

        mDbAccess = hasDatabaseConfig(mCmdLineArgs) ? createDatabaseAccess(mCmdLineArgs) : null;

        final int threads = mCmdLineArgs.hasOption(THREADS) ? Integer.parseInt(mCmdLineArgs.getOptionValue(THREADS)) : THREADS_DEFAULT;

        mExecutorService = Executors.newFixedThreadPool(threads);

        // load config
        mConfig = new PurpleConfig(mPurpleVersion.version(), mCmdLineArgs);

        // and common reference data
        mReferenceData = new ReferenceData(mCmdLineArgs, mConfig);

        mGermlineVariants = new GermlineVariants(mConfig, mReferenceData, mPurpleVersion.version());

        if(!mConfig.DriversOnly)
        {
            mSegmentation = new Segmentation(mConfig, mReferenceData);
            mCharts = new Charts(mConfig, mExecutorService, mReferenceData.RefGenVersion.is38());
        }
        else
        {
            mSegmentation = null;
            mCharts = null;
        }

        if(!mConfig.isValid() || !mReferenceData.isValid())
        {
            PPL_LOGGER.error("initialisation error, exiting");
            mExecutorService.shutdown();
            System.exit(1);
        }
    }

    public void run()
    {
        try
        {
            processSample(mConfig.ReferenceId, mConfig.TumorId);

        }
        finally
        {
            mExecutorService.shutdown();
        }

        PPL_LOGGER.info("Complete");
    }

    private SampleData loadSampleData(final String referenceId, final String tumorSample, final SampleDataFiles sampleDataFiles)
    {
        try
        {
            SampleData sampleData = null;

            if(!mConfig.DriversOnly)
            {
                // load amber and cobalt sample data
                final AmberData amberData = new AmberData(tumorSample, sampleDataFiles.AmberDirectory);

                final CobaltData cobaltData = new CobaltData(referenceId,
                        tumorSample,
                        sampleDataFiles.CobaltDirectory,
                        amberData.PatientGender,
                        mConfig.TumorOnlyMode);

                // load structural and somatic variants
                final StructuralVariantCache svCache = createStructuralVariantCache(tumorSample, sampleDataFiles);

                sampleData = new SampleData(referenceId, tumorSample, amberData, cobaltData, svCache);
            }
            else
            {
                sampleData = new SampleData(referenceId, tumorSample, null, null, null);
            }

            sampleData.loadSomatics(sampleDataFiles.SomaticVcfFile, mReferenceData, mConfig.TumorOnlyMode);

            return sampleData;
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed processing sample({}): {}", tumorSample, e.toString());
            return null;
        }
    }

    private void processSample(final String referenceId, final String tumorSample)
    {
        PPL_LOGGER.info("Processing sample(ref={} tumor={})", referenceId, tumorSample);

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
                    return;

                performFit(referenceId, tumorSample, sampleDataFiles, sampleData);
            }
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed processing sample({}): {}", tumorSample, e.toString());
            e.printStackTrace();
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
        final CobaltChromosomes cobaltChromosomes = cobaltData.CobaltChromosomes;

        if(cobaltGender.equals(amberGender))
        {
            PPL_LOGGER.info("Sample gender is {}", cobaltGender.toString().toLowerCase());
        }
        else
        {
            PPL_LOGGER.warn("COBALT gender {} does not match AMBER gender {}", cobaltGender, amberGender);
        }

        PPL_LOGGER.info("Applying segmentation");
        final List<ObservedRegion> observedRegions = mSegmentation.createSegments(sampleData.SvCache.variants(), amberData, cobaltData);

        PPL_LOGGER.info("Fitting purity");
        final FittedRegionFactory fittedRegionFactory =
                createFittedRegionFactory(amberData.AverageTumorDepth, cobaltChromosomes, mConfig.Fitting);

        final BestFit bestFit = fitPurity(sampleData, observedRegions, fittedRegionFactory, sampleData.SvCache.variants());

        final FittedPurity fittedPurity = bestFit.fit();

        final PurityAdjuster purityAdjuster =
                new PurityAdjusterAbnormalChromosome(fittedPurity.purity(), fittedPurity.normFactor(), cobaltChromosomes.chromosomes());

        final PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(
                mConfig.Fitting.MinDiploidTumorRatioCount,
                mConfig.Fitting.MinDiploidTumorRatioCountAtCentromere,
                amberData.AverageTumorDepth,
                fittedPurity.ploidy(),
                purityAdjuster,
                cobaltData.CobaltChromosomes);

        PPL_LOGGER.info("Calculating copy number");
        List<FittedRegion> fittedRegions =
                fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), observedRegions);
        copyNumberFactory.invoke(fittedRegions, sampleData.SvCache.variants());

        final int recoveredSVCount = recoverStructuralVariants(
                sampleData, sampleDataFiles, purityAdjuster, copyNumberFactory.copyNumbers());

        if(recoveredSVCount > 0)
        {
            PPL_LOGGER.info("Reapplying segmentation with {} recovered structural variants", recoveredSVCount);
            final List<ObservedRegion> recoveredObservedRegions =
                    mSegmentation.createSegments(sampleData.SvCache.variants(), amberData, cobaltData);

            PPL_LOGGER.info("Recalculating copy number");
            fittedRegions = fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), recoveredObservedRegions);
            copyNumberFactory.invoke(fittedRegions, sampleData.SvCache.variants());
        }

        final List<PurpleCopyNumber> copyNumbers = copyNumberFactory.copyNumbers();
        sampleData.SvCache.inferMissingVariant(copyNumbers);

        final List<FittedRegion> enrichedFittedRegions = updateRegionsWithCopyNumbers(fittedRegions, copyNumbers);

        final List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFactory.geneCopyNumbers(mReferenceData.GeneTransCache, copyNumbers);

        PPL_LOGGER.info("Generating QC Stats");
        final PurpleQC qcChecks = PurpleQCFactory.create(amberData.Contamination,
                bestFit,
                amberGender,
                cobaltGender,
                copyNumbers,
                geneCopyNumbers,
                cobaltChromosomes.germlineAberrations(),
                amberData.AverageTumorDepth);

        PPL_LOGGER.info("Modelling somatic peaks");
        final SomaticPeakStream somaticPeakStream = new SomaticPeakStream(mConfig);

        final List<PeakModel> somaticPeaks = somaticPeakStream.somaticPeakModel(
                purityAdjuster, copyNumbers, enrichedFittedRegions, sampleDataFiles.SomaticVcfFile);

        // at the moment the enriching of somatic variants is also contributing to the purity context,
        // so it cannot be done afterwards
        // if the read and write process were split then so could the fitting and enriching steps
        PPL_LOGGER.info("Enriching somatic variants");

        final SomaticStream somaticStream = new SomaticStream(
                mConfig, mReferenceData, somaticPeakStream.snpCount(), somaticPeakStream.indelCount(),
                somaticPeaks, sampleDataFiles.SomaticVcfFile);

        somaticStream.processAndWrite(purityAdjuster, copyNumbers, enrichedFittedRegions);

        final Set<String> reportedGenes = somaticStream.reportedGenes();

        if(!sampleDataFiles.GermlineVcfFile.isEmpty())
        {
            mGermlineVariants.processAndWrite(referenceId,
                    tumorSample, sampleDataFiles.GermlineVcfFile, purityAdjuster, copyNumbers, reportedGenes);
        }

        final PurityContext purityContext = ImmutablePurityContext.builder()
                .version(mPurpleVersion.version())
                .bestFit(bestFit.fit())
                .method(bestFit.method())
                .gender(cobaltGender)
                .score(bestFit.score())
                .polyClonalProportion(polyclonalProportion(copyNumbers))
                .wholeGenomeDuplication(wholeGenomeDuplication(copyNumbers))
                .microsatelliteIndelsPerMb(somaticStream.microsatelliteIndelsPerMb())
                .microsatelliteStatus(somaticStream.microsatelliteStatus())
                .tumorMutationalLoad(somaticStream.tumorMutationalLoad())
                .tumorMutationalLoadStatus(somaticStream.tumorMutationalLoadStatus())
                .tumorMutationalBurdenPerMb(somaticStream.tumorMutationalBurdenPerMb())
                .tumorMutationalBurdenStatus(somaticStream.tumorMutationalBurdenPerMbStatus())
                .svTumorMutationalBurden(sampleData.SvCache.passingBnd())
                .qc(qcChecks)
                .build();

        GermlineDeletionDrivers germlineDeletionDrivers = new GermlineDeletionDrivers(
                mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

        germlineDeletionDrivers.findDeletions(copyNumbers, fittedRegions);
        final List<GermlineDeletion> germlineDeletions = germlineDeletionDrivers.getDeletions();

        // write fit results
        PPL_LOGGER.info("Writing purple data to directory: {}", mConfig.OutputDir);
        mPurpleVersion.write(mConfig.OutputDir);
        PurityContextFile.write(mConfig.OutputDir, tumorSample, purityContext);
        FittedPurityRangeFile.write(mConfig.OutputDir, tumorSample, bestFit.allFits());
        PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample), copyNumbers);
        GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample), geneCopyNumbers);
        SegmentFile.write(SegmentFile.generateFilename(mConfig.OutputDir, tumorSample), fittedRegions);
        GermlineDeletion.write(GermlineDeletion.generateFilename(mConfig.OutputDir, tumorSample), germlineDeletions);

        if(mDbAccess != null)
        {
            PPL_LOGGER.info("Writing purple data to database: {}", mCmdLineArgs.getOptionValue(DB_URL));

            mDbAccess.writePurity(tumorSample, purityContext, qcChecks);
            mDbAccess.writeBestFitPerPurity(tumorSample, bestFit.bestFitPerPurity());
            mDbAccess.writeCopynumbers(tumorSample, copyNumbers);
            mDbAccess.writeGeneCopyNumbers(tumorSample, geneCopyNumbers);
            mDbAccess.writeGermlineDeletions(tumorSample, germlineDeletions);
        }

        if(mConfig.Charting.Enabled || mConfig.Charting.CircosBinary.isPresent())
        {
            PPL_LOGGER.info("Generating charts");

            mCharts.write(
                    referenceId, tumorSample, !sampleDataFiles.SomaticVcfFile.isEmpty(),
                    cobaltGender, copyNumbers, somaticStream.downsampledVariants(), sampleData.SvCache.variants(),
                    fittedRegions, Lists.newArrayList(amberData.ChromosomeBafs.values()));
        }

        if(mConfig.RunDrivers)
        {
            findDrivers(tumorSample, purityContext, copyNumbers, geneCopyNumbers, fittedRegions, somaticStream);
        }
    }

    private void findDrivers(final String tumorSample) throws IOException
    {
        final String purpleDataPath = mConfig.OutputDir;

        String somaticVcf = purpleDataPath + tumorSample + PURPLE_SOMATIC_VCF_SUFFIX;

        // the counts passed in here are only for down-sampling for charting, which is not relevant for drivers
        final SomaticStream somaticStream = new SomaticStream(
                mConfig, mReferenceData, 0, 0, null, somaticVcf);

        int somaticVariantCount = somaticStream.loadVariantsForDrivers(somaticVcf);
        PPL_LOGGER.info("loaded {} somatic variants", somaticVariantCount);

        final String germlineVcf = purpleDataPath + tumorSample + PURPLE_GERMLINE_VCF_SUFFIX;

        if(Files.exists(Paths.get(germlineVcf)))
        {
            mGermlineVariants.loadReportableVariants(germlineVcf);
            PPL_LOGGER.info("loaded {} reportable germline variants", mGermlineVariants.reportableVariants().size());
        }

        final PurityContext purityContext = PurityContextFile.read(purpleDataPath, tumorSample);

        final List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                PurpleCopyNumberFile.generateFilenameForReading(purpleDataPath, tumorSample));

        final List<FittedRegion> fittedRegions = SegmentFile.read(SegmentFile.generateFilename(purpleDataPath, tumorSample)).stream()
                .filter(x -> x.germlineStatus() == HET_DELETION || x.germlineStatus() == HOM_DELETION)
                .collect(Collectors.toList());

        final List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(
                GeneCopyNumberFile.generateFilenameForReading(purpleDataPath, tumorSample));

        findDrivers(tumorSample, purityContext, copyNumbers, geneCopyNumbers, fittedRegions, somaticStream);
    }

    private void findDrivers(
            final String tumorSample, final PurityContext purityContext, final List<PurpleCopyNumber> copyNumbers,
            final List<GeneCopyNumber> geneCopyNumbers, final List<FittedRegion> fittedRegions,
            final SomaticStream somaticStream) throws IOException
    {
        final List<DriverCatalog> somaticDriverCatalog = Lists.newArrayList();
        final List<DriverCatalog> germlineDriverCatalog = Lists.newArrayList();

        PPL_LOGGER.info("Generating driver catalog");

        final Map<String,List<GeneCopyNumber>> geneCopyNumberMap = listToMap(geneCopyNumbers);

        somaticDriverCatalog.addAll(somaticStream.drivers(geneCopyNumberMap));

        final CNADrivers cnaDrivers = new CNADrivers(purityContext.qc().status(), mReferenceData.DriverGenes);

        somaticDriverCatalog.addAll(cnaDrivers.deletions(geneCopyNumbers));
        somaticDriverCatalog.addAll(cnaDrivers.amplifications(purityContext.bestFit().ploidy(), geneCopyNumbers));

        final GermlineDrivers germlineDrivers = new GermlineDrivers(mReferenceData.DriverGenes.driverGenes());
        germlineDriverCatalog.addAll(germlineDrivers.findDrivers(mGermlineVariants.reportableVariants(), geneCopyNumberMap));

        GermlineDeletionDrivers germlineDeletionDrivers = new GermlineDeletionDrivers(
                mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

        germlineDeletionDrivers.findDeletions(copyNumbers, fittedRegions);

        germlineDriverCatalog.addAll(germlineDeletionDrivers.getDrivers());

        DriverCatalogFile.write(DriverCatalogFile.generateSomaticFilename(mConfig.OutputDir, tumorSample), somaticDriverCatalog);
        DriverCatalogFile.write(DriverCatalogFile.generateGermlineFilename(mConfig.OutputDir, tumorSample), germlineDriverCatalog);

        if(mDbAccess != null)
            mDbAccess.writePurpleDriverCatalog(tumorSample, somaticDriverCatalog, germlineDriverCatalog);
    }

    private int recoverStructuralVariants(final SampleData sampleData, final SampleDataFiles sampleDataFiles,
            final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        if(sampleDataFiles.RecoveredSvVcfFile.isEmpty())
        {
            return 0;
        }

        PPL_LOGGER.info("Loading recovery candidates from {}", sampleDataFiles.RecoveredSvVcfFile);

        try(final RecoverStructuralVariants recovery = new RecoverStructuralVariants(mConfig.Misc,
                purityAdjuster,
                sampleDataFiles.RecoveredSvVcfFile,
                copyNumbers))
        {
            final Collection<VariantContext> recoveredVariants = recovery.recoverVariants(sampleData.SvCache.variants());

            if(!recoveredVariants.isEmpty())
            {
                recoveredVariants.forEach(x -> sampleData.SvCache.addVariant(x));
            }

            return recoveredVariants.size();
        }
    }

    private BestFit fitPurity(final SampleData sampleData, final List<ObservedRegion> observedRegions,
            final FittedRegionFactory fittedRegionFactory, final List<StructuralVariant> structuralVariants)
            throws ExecutionException, InterruptedException
    {
        final FittingConfig fittingConfig = mConfig.Fitting;
        final SomaticFitConfig somaticFitConfig = mConfig.SomaticFitting;

        final List<FittedPurity> fitCandidates = Lists.newArrayList();

        final CobaltChromosomes cobaltChromosomes = sampleData.Cobalt.CobaltChromosomes;
        final List<SomaticVariant> snpSomatics = sampleData.FittingSomaticVariants;

        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(mExecutorService,
                cobaltChromosomes,
                fittingConfig.MinPurity,
                fittingConfig.MaxPurity,
                fittingConfig.PurityIncrement,
                fittingConfig.MinPloidy,
                fittingConfig.MaxPloidy,
                somaticFitConfig.SomaticPenaltyWeight,
                mConfig.TumorOnlyMode,
                fittedRegionFactory,
                observedRegions,
                snpSomatics);

        fitCandidates.addAll(fittedPurityFactory.all());

        final BestFitFactory bestFitFactory = new BestFitFactory(mConfig,
                sampleData.Amber.minSomaticTotalReadCount(),
                sampleData.Amber.maxSomaticTotalReadCount(),
                fitCandidates,
                snpSomatics,
                structuralVariants,
                observedRegions);

        return bestFitFactory.bestFit();
    }

    private FittedRegionFactory createFittedRegionFactory(final int averageTumorDepth, final CobaltChromosomes cobaltChromosomes,
            final FittingConfig fitScoreConfig)
    {
        return new FittedRegionFactoryV2(cobaltChromosomes,
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

        PPL_LOGGER.info("Loading structural variants from {}", sampleDataFiles.SvVcfFile);

        final String outputVcf = mConfig.OutputDir + File.separator + tumorSample + ".purple.sv.vcf.gz";

        return new StructuralVariantCache(mPurpleVersion.version(), sampleDataFiles.SvVcfFile, outputVcf, mReferenceData);
    }

    public static void main(final String... args) throws IOException
    {
        final Options options = createOptions();

        try
        {
            PurpleApplication purpleApplication = new PurpleApplication(options, args);
            purpleApplication.run();
        } catch(ParseException e)
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
