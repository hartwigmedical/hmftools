package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

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
import com.hartwig.hmftools.common.amber.BaseDepthEvidence;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.BaseDepthFilter;
import com.hartwig.hmftools.common.amber.NormalHeterozygousFilter;
import com.hartwig.hmftools.common.amber.NormalHomozygousFilter;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAFEvidence;
import com.hartwig.hmftools.common.amber.TumorContamination;
import com.hartwig.hmftools.common.amber.TumorContaminationEvidence;
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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(AmberApplication.class);

    private final AmberConfig config;
    private final ExecutorService executorService;
    private final Predicate<BaseDepth> snpCheckFilter;
    private final Predicate<BaseDepth> homozygousFilter;
    private final Predicate<BaseDepth> heterozygousFilter;
    private final AmberPersistence persistence;
    private final VersionInfo versionInfo;
    private final ListMultimap<Chromosome, AmberSite> sites;

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException {
        final Options options = AmberConfig.createOptions();
        try (final AmberApplication application = new AmberApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("AmberApplication", options);
            System.exit(1);
        }
    }

    private AmberApplication(final Options options, final String... args) throws IOException, ParseException {
        versionInfo = new VersionInfo("amber.version");
        LOGGER.info("AMBER version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        config = AmberConfig.createConfig(cmd);

        final Predicate<BaseDepth> isValidFilter = BaseDepth::isValid;
        homozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        heterozygousFilter = new NormalHeterozygousFilter(config.minHetAfPercent(), config.maxHetAfPercent()).and(isValidFilter);
        persistence = new AmberPersistence(config);

        final File outputDir = new File(config.outputDirectory());
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + config.outputDirectory());
        }

        if (!new File(config.bafLociPath()).exists()) {
            throw new IOException("Unable to locate vcf file " + config.bafLociPath());
        }

        if (!new File(config.tumorBamPath()).exists()) {
            throw new IOException("Unable to locate tumor bam file " + config.tumorBamPath());
        }

        if (!config.tumorOnly()) {
            for (String referenceBam : config.referenceBamPath()) {
                if (!new File(referenceBam).exists()) {
                    throw new ParseException("Unable to locate reference bam " + referenceBam);
                }
            }
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        executorService = Executors.newFixedThreadPool(config.threadCount(), namedThreadFactory);

        LOGGER.info("Loading vcf file {}", config.bafLociPath());
        sites = AmberSiteFactory.sites(config.bafLociPath());
        snpCheckFilter = new SnpCheckFilter(sites);
    }

    private void run() throws InterruptedException, ExecutionException, IOException {
        if (config.tumorOnly()) {
            runTumorOnly();
        } else {
            runNormalMode();
        }
    }

    private void runNormalMode() throws InterruptedException, ExecutionException, IOException {
        final SamReaderFactory readerFactory = readerFactory(config);
        final AmberHetNormalEvidence hetNormalEvidence = new AmberHetNormalEvidence();

        // Primary Reference Data
        final ListMultimap<Chromosome, BaseDepth> unfilteredNormal = normalDepth(readerFactory, config.referenceBamPath().get(0), sites);
        final Predicate<BaseDepth> depthFilter = new BaseDepthFilter(config.minDepthPercent(), config.maxDepthPercent(), unfilteredNormal);
        final ListMultimap<Chromosome, BaseDepth> snpCheck = filterEntries(unfilteredNormal, snpCheckFilter);
        final ListMultimap<Chromosome, BaseDepth> homNormal = filterEntries(unfilteredNormal, depthFilter.and(homozygousFilter));
        final ListMultimap<Chromosome, BaseDepth> hetNormal = filterEntries(unfilteredNormal, depthFilter.and(heterozygousFilter));
        hetNormalEvidence.add(config.primaryReference(), hetNormal.values());

        // Additional Reference Data
        for (int i = 1; i < config.reference().size(); i++) {
            final String sample = config.reference().get(i);
            final String sampleBam = config.referenceBamPath().get(i);
            final Collection<BaseDepth> additional = normalDepth(readerFactory, sampleBam, hetNormalEvidence.intersection()).values();
            final Predicate<BaseDepth> filter = new BaseDepthFilter(config.minDepthPercent(), config.maxDepthPercent(), additional);
            final Collection<BaseDepth> additionalHetNormal = additional.stream().filter(filter.and(heterozygousFilter)).collect(toList());
            hetNormalEvidence.add(sample, additionalHetNormal);
        }

        final Predicate<BaseDepth> intersectionFilter = hetNormalEvidence.intersectionFilter();
        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, filterEntries(hetNormal, intersectionFilter));
        final List<TumorBAF> tumorBAFList = tumorBAFMap.values().stream().sorted().collect(toList());
        final List<AmberBAF> amberBAFList = tumorBAFList.stream().map(AmberBAF::create).filter(AmberApplication::isValid).collect(toList());

        final ListMultimap<Chromosome, TumorContamination> tumorContamination = contamination(readerFactory, homNormal);
        final List<TumorContamination> contaminationList = Lists.newArrayList(tumorContamination.values());

        persistence.persisQC(amberBAFList, contaminationList);
        persistence.persistVersionInfo(versionInfo);
        persistence.persistBafVcf(tumorBAFList, hetNormalEvidence);
        persistence.persistContamination(contaminationList);
        persistence.persistSnpCheck(snpCheck);
        persistence.persistBAF(amberBAFList);
    }

    private void runTumorOnly() throws InterruptedException, ExecutionException, IOException {
        final SamReaderFactory readerFactory = readerFactory(config);

        final ListMultimap<Chromosome, BaseDepth> allNormal = emptyNormalHetSites(sites);
        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, allNormal);

        final List<TumorBAF> tumorBAFList = tumorBAFMap.values()
                .stream()
                .filter(x -> x.tumorRefSupport() >= config.tumorOnlyMinSupport())
                .filter(x -> x.tumorAltSupport() >= config.tumorOnlyMinSupport())
                .filter(x -> isFinite(x.refFrequency()) && Doubles.greaterOrEqual(x.refFrequency(), config.tumorOnlyMinVaf()))
                .filter(x -> isFinite(x.altFrequency()) && Doubles.greaterOrEqual(x.altFrequency(), config.tumorOnlyMinVaf()))
                .sorted()
                .collect(toList());
        final List<AmberBAF> amberBAFList =
                tumorBAFList.stream().map(AmberBAF::create).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        persistence.persisQC(amberBAFList, Lists.newArrayList());
        persistence.persistVersionInfo(versionInfo);
        persistence.persistBafVcf(tumorBAFList, new AmberHetNormalEvidence());
        persistence.persistBAF(amberBAFList);
    }

    @NotNull
    private ListMultimap<Chromosome, BaseDepth> normalDepth(final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome, AmberSite> bedRegionsSortedSet) throws InterruptedException, ExecutionException {

        final int partitionSize = Math.max(config.minPartition(), bedRegionsSortedSet.size() / config.threadCount());

        LOGGER.info("Processing {} potential sites in reference bam {}", bedRegionsSortedSet.values().size(), bamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<BaseDepthEvidence>> futures = Lists.newArrayList();
        for (final Chromosome contig : bedRegionsSortedSet.keySet()) {
            for (final List<AmberSite> inner : Lists.partition(Lists.newArrayList(bedRegionsSortedSet.get(contig)), partitionSize)) {
                final BaseDepthEvidence evidence = new BaseDepthEvidence(config.typicalReadDepth(),
                        config.minMappingQuality(),
                        config.minBaseQuality(),
                        inner.get(0).chromosome(),
                        bamPath,
                        readerFactory,
                        inner);
                futures.add(executorService.submit(completion.task(evidence)));
            }
        }

        final ListMultimap<Chromosome, BaseDepth> normalEvidence = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));
        return normalEvidence;
    }

    @NotNull
    private ListMultimap<Chromosome, BaseDepth> emptyNormalHetSites(@NotNull final ListMultimap<Chromosome, AmberSite> sites) {
        final ListMultimap<Chromosome, BaseDepth> result = ArrayListMultimap.create();
        for (Chromosome chromosome : sites.keySet()) {
            result.putAll(chromosome, sites.get(chromosome).stream().map(BaseDepthFactory::fromAmberSite).collect(toList()));
        }

        return result;
    }

    @NotNull
    private ListMultimap<Chromosome, TumorBAF> tumorBAF(@NotNull final SamReaderFactory readerFactory,
            @NotNull final ListMultimap<Chromosome, BaseDepth> normalHetSites) throws ExecutionException, InterruptedException {
        final int partitionSize = Math.max(config.minPartition(), normalHetSites.values().size() / config.threadCount());

        LOGGER.info("Processing {} heterozygous sites in tumor bam {}", normalHetSites.values().size(), config.tumorBamPath());
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorBAFEvidence>> futures = Lists.newArrayList();
        for (final Chromosome chromosome : normalHetSites.keySet()) {
            for (final List<BaseDepth> chromosomeBafPoints : Lists.partition(normalHetSites.get(chromosome), partitionSize)) {
                if (!chromosomeBafPoints.isEmpty()) {
                    final String contig = chromosomeBafPoints.get(0).chromosome();
                    final TumorBAFEvidence evidence = new TumorBAFEvidence(config.typicalReadDepth(),
                            config.minMappingQuality(),
                            config.minBaseQuality(),
                            contig,
                            config.tumorBamPath(),
                            readerFactory,
                            chromosomeBafPoints);
                    futures.add(executorService.submit(completion.task(evidence)));
                }
            }
        }

        final ListMultimap<Chromosome, TumorBAF> result = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    @NotNull
    private ListMultimap<Chromosome, TumorContamination> contamination(@NotNull final SamReaderFactory readerFactory,
            @NotNull final ListMultimap<Chromosome, BaseDepth> normalHomSites) throws ExecutionException, InterruptedException {
        final int partitionSize = Math.max(config.minPartition(), normalHomSites.values().size() / config.threadCount());

        LOGGER.info("Processing {} homozygous sites in tumor bam {} for contamination", normalHomSites.size(), config.tumorBamPath());
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorContaminationEvidence>> futures = Lists.newArrayList();
        for (final Chromosome chromosome : normalHomSites.keySet()) {
            for (final List<BaseDepth> chromosomeBafPoints : Lists.partition(normalHomSites.get(chromosome), partitionSize)) {
                if (!chromosomeBafPoints.isEmpty()) {
                    final String contig = chromosomeBafPoints.get(0).chromosome();
                    final TumorContaminationEvidence evidence = new TumorContaminationEvidence(config.typicalReadDepth(),
                            config.minMappingQuality(),
                            config.minBaseQuality(),
                            contig,
                            config.tumorBamPath(),
                            readerFactory,
                            chromosomeBafPoints);
                    futures.add(executorService.submit(completion.task(evidence)));
                }
            }
        }

        final ListMultimap<Chromosome, TumorContamination> result = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static boolean isValid(@NotNull final AmberBAF baf) {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
    }

    @NotNull
    private static <T> List<T> getFuture(@NotNull final List<Future<T>> futures) throws ExecutionException, InterruptedException {
        final List<T> result = Lists.newArrayList();
        for (Future<T> chromosomeBAFEvidenceFuture : futures) {
            result.add(chromosomeBAFEvidenceFuture.get());
        }
        return result;
    }

    @NotNull
    private static SamReaderFactory readerFactory(@NotNull final AmberConfig config) {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.validationStringency());
        if (!config.refGenomePath().isEmpty()) {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.refGenomePath())));
        }
        return readerFactory;
    }

    @Override
    public void close() {
        executorService.shutdown();
        LOGGER.info("Complete");
    }
}
