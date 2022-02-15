package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.collection.Multimaps.filterEntries;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.NormalHeterozygousFilter;
import com.hartwig.hmftools.common.amber.NormalHomozygousFilter;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable
{
    // add the AmberConfig options
    @ParametersDelegate
    private final AmberConfig mConfig = new AmberConfig();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    private ExecutorService mExecutorService;
    private Predicate<BaseDepth> mSnpCheckFilter;
    private Predicate<BaseDepth> mHomozygousFilter;
    private Predicate<BaseDepth> mHeterozygousFilter;
    private AmberPersistence mPersistence;
    private VersionInfo mVersionInfo;
    private ListMultimap<Chromosome,AmberSite> mChromosomeSites;

    public int run() throws IOException, InterruptedException, ExecutionException
    {
        mLoggingOptions.setLogLevel();

        mVersionInfo = new VersionInfo("amber.version");
        AMB_LOGGER.info("AMBER version: {}", mVersionInfo.version());

        final Predicate<BaseDepth> isValidFilter = BaseDepth::isValid;
        mHomozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        mHeterozygousFilter = new NormalHeterozygousFilter(mConfig.MinHetAfPercent, mConfig.MaxHetAfPercent).and(isValidFilter);
        mPersistence = new AmberPersistence(mConfig);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        AMB_LOGGER.info("Loading vcf file {}", mConfig.BafLociPath);
        mChromosomeSites = AmberSiteFactory.sites(mConfig.BafLociPath);
        mSnpCheckFilter = new SnpCheckFilter(mChromosomeSites);

        if(!mConfig.isValid())
        {
            AMB_LOGGER.error(" invalid config, exiting");
            return 1;
        }

        if(mConfig.isTumorOnly())
        {
            runTumorOnly();
        }
        else
        {
            runNormalMode();
        }
        return 0;
    }

    private void runNormalMode() throws InterruptedException, ExecutionException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);
        final AmberHetNormalEvidence hetNormalEvidence = new AmberHetNormalEvidence();

        // Primary Reference Data
        final ListMultimap<Chromosome, BaseDepth> unfilteredNormal = normalDepth(readerFactory, mConfig.ReferenceBamPath.get(0), mChromosomeSites);
        final Predicate<BaseDepth> depthFilter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, unfilteredNormal);
        final ListMultimap<Chromosome, BaseDepth> snpCheck = filterEntries(unfilteredNormal, mSnpCheckFilter);
        final ListMultimap<Chromosome, BaseDepth> homNormal = filterEntries(unfilteredNormal, depthFilter.and(mHomozygousFilter));
        final ListMultimap<Chromosome, BaseDepth> hetNormal = filterEntries(unfilteredNormal, depthFilter.and(mHeterozygousFilter));
        hetNormalEvidence.add(mConfig.primaryReference(), hetNormal.values());

        // Additional Reference Data
        for(int i = 1; i < mConfig.ReferenceIds.size(); i++)
        {
            final String sample = mConfig.ReferenceIds.get(i);
            final String sampleBam = mConfig.ReferenceBamPath.get(i);
            final Collection<BaseDepth> additional = normalDepth(readerFactory, sampleBam, hetNormalEvidence.intersection()).values();
            final Predicate<BaseDepth> filter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, additional);
            final Collection<BaseDepth> additionalHetNormal = additional.stream().filter(filter.and(mHeterozygousFilter)).collect(toList());
            hetNormalEvidence.add(sample, additionalHetNormal);
        }

        final Predicate<BaseDepth> intersectionFilter = hetNormalEvidence.intersectionFilter();
        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, filterEntries(hetNormal, intersectionFilter));
        final List<TumorBAF> tumorBAFList = tumorBAFMap.values().stream().sorted().collect(toList());
        final List<AmberBAF> amberBAFList = tumorBAFList.stream().map(AmberBAF::create).filter(AmberApplication::isValid).collect(toList());

        final ListMultimap<Chromosome, TumorContamination> tumorContamination = contamination(readerFactory, homNormal);
        final List<TumorContamination> contaminationList = new ArrayList<>(tumorContamination.values());

        RegionOfHomozygosityFinder rohFinder = new RegionOfHomozygosityFinder(mConfig.refGenomeVersion, mConfig.MinDepthPercent, mConfig.MaxDepthPercent);
        final List<RegionOfHomozygosity> regionsOfHomozygosity = rohFinder.findRegions(unfilteredNormal);

        double consanguinityProportion = ConsanguinityAnalyser.calcConsanguinityProportion(regionsOfHomozygosity);
        Chromosome uniparentalDisomy = ConsanguinityAnalyser.findUniparentalDisomy(regionsOfHomozygosity);

        mPersistence.persistQC(amberBAFList, contaminationList, consanguinityProportion, uniparentalDisomy);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBafVcf(tumorBAFList, hetNormalEvidence);
        mPersistence.persistContamination(contaminationList);
        mPersistence.persistSnpCheck(snpCheck);
        //mPersistence.persistPrimaryRefUnfiltered(unfilteredNormal);
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(regionsOfHomozygosity);
    }

    private void runTumorOnly() throws InterruptedException, ExecutionException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);

        final ListMultimap<Chromosome, BaseDepth> allNormal = emptyNormalHetSites(mChromosomeSites);
        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, allNormal);

        final List<TumorBAF> tumorBAFList = tumorBAFMap.values()
                .stream()
                .filter(x -> x.tumorRefSupport() >= mConfig.TumorOnlyMinSupport)
                .filter(x -> x.tumorAltSupport() >= mConfig.TumorOnlyMinSupport)
                .filter(x -> isFinite(x.refFrequency()) && Doubles.greaterOrEqual(x.refFrequency(), mConfig.TumorOnlyMinVaf))
                .filter(x -> isFinite(x.altFrequency()) && Doubles.greaterOrEqual(x.altFrequency(), mConfig.TumorOnlyMinVaf))
                .sorted()
                .collect(toList());
        final List<AmberBAF> amberBAFList =
                tumorBAFList.stream().map(AmberBAF::create).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        mPersistence.persistQC(amberBAFList, new ArrayList<>(), 0.0, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBafVcf(tumorBAFList, new AmberHetNormalEvidence());
        mPersistence.persistBAF(amberBAFList);
    }

    private ListMultimap<Chromosome, BaseDepth> normalDepth(
            final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome, AmberSite> bedRegionsSortedSet) throws InterruptedException, ExecutionException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), bedRegionsSortedSet.size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} potential sites in reference bam {}", bedRegionsSortedSet.values().size(), bamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<BaseDepthEvidence>> futures = new ArrayList<>();
        for(final Chromosome contig : bedRegionsSortedSet.keySet())
        {
            for(final List<AmberSite> inner : Lists.partition(new ArrayList<>(bedRegionsSortedSet.get(contig)), partitionSize))
            {
                final BaseDepthEvidence evidence = new BaseDepthEvidence(mConfig.typicalReadDepth(),
                        mConfig.MinMappingQuality, mConfig.MinBaseQuality,
                        inner.get(0).chromosome(), bamPath, readerFactory, inner);
                
                futures.add(mExecutorService.submit(completion.task(evidence)));
            }
        }

        final ListMultimap<Chromosome, BaseDepth> normalEvidence = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));
        return normalEvidence;
    }

    private ListMultimap<Chromosome, BaseDepth> emptyNormalHetSites(final ListMultimap<Chromosome, AmberSite> sites)
    {
        final ListMultimap<Chromosome, BaseDepth> result = ArrayListMultimap.create();
        for(Chromosome chromosome : sites.keySet())
        {
            result.putAll(chromosome, sites.get(chromosome).stream().map(BaseDepthFactory::fromAmberSite).collect(toList()));
        }

        return result;
    }

    private ListMultimap<Chromosome, TumorBAF> tumorBAF(final SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, BaseDepth> normalHetSites) throws ExecutionException, InterruptedException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), normalHetSites.values().size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} heterozygous sites in tumor bam {}", normalHetSites.values().size(), mConfig.TumorBamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorBAFEvidence>> futures = new ArrayList<>();
        for(final Chromosome chromosome : normalHetSites.keySet())
        {
            for(final List<BaseDepth> chromosomeBafPoints : Lists.partition(normalHetSites.get(chromosome), partitionSize))
            {
                if(!chromosomeBafPoints.isEmpty())
                {
                    final String contig = chromosomeBafPoints.get(0).chromosome();
                    final TumorBAFEvidence evidence = new TumorBAFEvidence(mConfig.typicalReadDepth(),
                            mConfig.MinMappingQuality,
                            mConfig.MinBaseQuality,
                            contig,
                            mConfig.TumorBamPath,
                            readerFactory,
                            chromosomeBafPoints);
                    futures.add(mExecutorService.submit(completion.task(evidence)));
                }
            }
        }

        final ListMultimap<Chromosome, TumorBAF> result = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    private ListMultimap<Chromosome, TumorContamination> contamination(
            final SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, BaseDepth> normalHomSites) throws ExecutionException, InterruptedException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), normalHomSites.values().size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} homozygous sites in tumor bam {} for contamination", normalHomSites.size(), mConfig.TumorBamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorContaminationEvidence>> futures = new ArrayList<>();
        for(final Chromosome chromosome : normalHomSites.keySet())
        {
            for(final List<BaseDepth> chromosomeBafPoints : Lists.partition(normalHomSites.get(chromosome), partitionSize))
            {
                if(!chromosomeBafPoints.isEmpty())
                {
                    final String contig = chromosomeBafPoints.get(0).chromosome();
                    final TumorContaminationEvidence evidence = new TumorContaminationEvidence(mConfig.typicalReadDepth(),
                            mConfig.MinMappingQuality,
                            mConfig.MinBaseQuality,
                            contig,
                            mConfig.TumorBamPath,
                            readerFactory,
                            chromosomeBafPoints);
                    futures.add(mExecutorService.submit(completion.task(evidence)));
                }
            }
        }

        final ListMultimap<Chromosome, TumorContamination> result = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    private static boolean isValid(final AmberBAF baf)
    {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
    }

    private static <T> List<T> getFuture(final List<Future<T>> futures) throws ExecutionException, InterruptedException
    {
        final List<T> result = new ArrayList<>();
        for(Future<T> chromosomeBAFEvidenceFuture : futures)
        {
            result.add(chromosomeBAFEvidenceFuture.get());
        }
        return result;
    }

    private static SamReaderFactory readerFactory(final AmberConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.Stringency);
        if(config.RefGenomePath != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomePath)));
        }
        return readerFactory;
    }

    @Override
    public void close()
    {
        mExecutorService.shutdown();
        AMB_LOGGER.info("Complete");
    }

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException
    {
        AmberApplication amberApp = new AmberApplication();
        JCommander commander = JCommander.newBuilder()
                .addObject(amberApp)
                .build();

        // use unix style formatter
        commander.setUsageFormatter(new UnixStyleUsageFormatter(commander));
        // help message show in order parameters are declared
        commander.setParameterDescriptionComparator(new DeclaredOrderParameterComparator(AmberApplication.class));

        try
        {
            commander.parse(args);
        }
        catch (com.beust.jcommander.ParameterException e)
        {
            System.out.println("Unable to parse args: " + e.getMessage());
            commander.usage();
            System.exit(1);
        }

        System.exit(amberApp.run());
    }

}
