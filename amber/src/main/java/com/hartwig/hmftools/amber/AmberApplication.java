package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.collection.Multimaps.filterEntries;

import java.io.File;
import java.io.IOException;
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
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable
{
    private final AmberConfig mConfig;
    private final ExecutorService mExecutorService;
    private final Predicate<BaseDepth> mSnpCheckFilter;
    private final Predicate<BaseDepth> mHomozygousFilter;
    private final Predicate<BaseDepth> mHeterozygousFilter;
    private final AmberPersistence mPersistence;
    private final VersionInfo mVersionInfo;
    private final ListMultimap<Chromosome,AmberSite> mChromosomeSites;

    private AmberApplication(final Options options, final String... args) throws IOException, ParseException
    {
        mVersionInfo = new VersionInfo("amber.version");
        AMB_LOGGER.info("AMBER version: {}", mVersionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        mConfig = new AmberConfig(cmd);

        final Predicate<BaseDepth> isValidFilter = BaseDepth::isValid;
        mHomozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        mHeterozygousFilter = new NormalHeterozygousFilter(mConfig.MinHetAfPercent, mConfig.MaxHetAfPercent).and(isValidFilter);
        mPersistence = new AmberPersistence(mConfig);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        AMB_LOGGER.info("Loading vcf file {}", mConfig.BafLociPath);
        mChromosomeSites = AmberSiteFactory.sites(mConfig.BafLociPath);
        mSnpCheckFilter = new SnpCheckFilter(mChromosomeSites);
    }

    private void run() throws InterruptedException, ExecutionException, IOException
    {
        if(!mConfig.isValid())
        {
            AMB_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }
        
        if(mConfig.TumorOnly)
        {
            runTumorOnly();
        }
        else
        {
            runNormalMode();
        }
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
        final List<TumorContamination> contaminationList = Lists.newArrayList(tumorContamination.values());

        mPersistence.persistQC(amberBAFList, contaminationList);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBafVcf(tumorBAFList, hetNormalEvidence);
        mPersistence.persistContamination(contaminationList);
        mPersistence.persistSnpCheck(snpCheck);
        mPersistence.persistBAF(amberBAFList);
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

        mPersistence.persistQC(amberBAFList, Lists.newArrayList());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBafVcf(tumorBAFList, new AmberHetNormalEvidence());
        mPersistence.persistBAF(amberBAFList);
    }

    @NotNull
    private ListMultimap<Chromosome, BaseDepth> normalDepth(final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome, AmberSite> bedRegionsSortedSet) throws InterruptedException, ExecutionException
    {

        final int partitionSize = Math.max(mConfig.minPartition(), bedRegionsSortedSet.size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} potential sites in reference bam {}", bedRegionsSortedSet.values().size(), bamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<BaseDepthEvidence>> futures = Lists.newArrayList();
        for(final Chromosome contig : bedRegionsSortedSet.keySet())
        {
            for(final List<AmberSite> inner : Lists.partition(Lists.newArrayList(bedRegionsSortedSet.get(contig)), partitionSize))
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

    @NotNull
    private ListMultimap<Chromosome, BaseDepth> emptyNormalHetSites(@NotNull final ListMultimap<Chromosome, AmberSite> sites)
    {
        final ListMultimap<Chromosome, BaseDepth> result = ArrayListMultimap.create();
        for(Chromosome chromosome : sites.keySet())
        {
            result.putAll(chromosome, sites.get(chromosome).stream().map(BaseDepthFactory::fromAmberSite).collect(toList()));
        }

        return result;
    }

    @NotNull
    private ListMultimap<Chromosome, TumorBAF> tumorBAF(@NotNull final SamReaderFactory readerFactory,
            @NotNull final ListMultimap<Chromosome, BaseDepth> normalHetSites) throws ExecutionException, InterruptedException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), normalHetSites.values().size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} heterozygous sites in tumor bam {}", normalHetSites.values().size(), mConfig.TumorBamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorBAFEvidence>> futures = Lists.newArrayList();
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

    @NotNull
    private ListMultimap<Chromosome, TumorContamination> contamination(@NotNull final SamReaderFactory readerFactory,
            @NotNull final ListMultimap<Chromosome, BaseDepth> normalHomSites) throws ExecutionException, InterruptedException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), normalHomSites.values().size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} homozygous sites in tumor bam {} for contamination", normalHomSites.size(), mConfig.TumorBamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorContaminationEvidence>> futures = Lists.newArrayList();
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

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static boolean isValid(@NotNull final AmberBAF baf)
    {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
    }

    @NotNull
    private static <T> List<T> getFuture(@NotNull final List<Future<T>> futures) throws ExecutionException, InterruptedException
    {
        final List<T> result = Lists.newArrayList();
        for(Future<T> chromosomeBAFEvidenceFuture : futures)
        {
            result.add(chromosomeBAFEvidenceFuture.get());
        }
        return result;
    }

    @NotNull
    private static SamReaderFactory readerFactory(@NotNull final AmberConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.Stringency);
        if(!config.RefGenomePath.isEmpty())
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
        final Options options = AmberConfig.createOptions();

        try
        {
            AmberApplication application = new AmberApplication(options, args);
            application.run();
        }
        catch(ParseException e)
        {
            AMB_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("AmberApplication", options);
            System.exit(1);
        }
    }

}
