package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

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
    private AmberPersistence mPersistence;
    private VersionInfo mVersionInfo;
    private ListMultimap<Chromosome,AmberSite> mChromosomeSites;

    public int run() throws IOException, InterruptedException, ExecutionException
    {
        mLoggingOptions.setLogLevel();

        mVersionInfo = new VersionInfo("amber.version");
        AMB_LOGGER.info("AMBER version: {}", mVersionInfo.version());

        mPersistence = new AmberPersistence(mConfig);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        AMB_LOGGER.info("Loading vcf file {}", mConfig.BafLociPath);
        mChromosomeSites = AmberSiteFactory.sites(mConfig.BafLociPath);

        if(!mConfig.isValid())
        {
            AMB_LOGGER.error(" invalid config, exiting");
            return 1;
        }

        if(mConfig.isTumorOnly())
        {
            runTumorOnly();
        }
        else if (mConfig.isGermlineOnly())
        {
            runGermlineOnly();
        }
        else
        {
            runNormalMode();
        }
        return 0;
    }

    private void runGermlineOnly() throws InterruptedException, ExecutionException, IOException
    {
        GermlineProcessor germlineProcessor = new GermlineProcessor(mConfig, readerFactory(mConfig), mExecutorService, mChromosomeSites);

        final List<AmberBAF> amberBAFList = germlineProcessor.getHeterozygousLoci().values().stream().map(AmberBAF::create)
                .filter(AmberApplication::isValid).sorted().collect(toList());

        mPersistence.persistQC(Collections.emptyList(), germlineProcessor.getConsanguinityProportion(), germlineProcessor.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistSnpCheck(germlineProcessor.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germlineProcessor.getRegionsOfHomozygosity());
    }

    private void runNormalMode() throws InterruptedException, ExecutionException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);

        GermlineProcessor germlineProcessor = new GermlineProcessor(mConfig, readerFactory(mConfig), mExecutorService, mChromosomeSites);

        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, germlineProcessor.getHeterozygousLoci());
        final List<TumorBAF> tumorBAFList = tumorBAFMap.values().stream().sorted().collect(toList());
        final List<AmberBAF> amberBAFList = tumorBAFList.stream().map(AmberBAF::create).filter(AmberApplication::isValid).collect(toList());

        final ListMultimap<Chromosome, TumorContamination> tumorContamination = contamination(readerFactory,
                germlineProcessor.getHomozygousLoci());
        final List<TumorContamination> contaminationList = new ArrayList<>(tumorContamination.values());

        mPersistence.persistQC(contaminationList, germlineProcessor.getConsanguinityProportion(), germlineProcessor.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistContamination(contaminationList);
        mPersistence.persistSnpCheck(germlineProcessor.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germlineProcessor.getRegionsOfHomozygosity());
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

        mPersistence.persistQC(Collections.emptyList(), 0.0, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBAF(amberBAFList);
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
        AmberUtils.getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

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
        AmberUtils.getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    private static boolean isValid(final AmberBAF baf)
    {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
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
