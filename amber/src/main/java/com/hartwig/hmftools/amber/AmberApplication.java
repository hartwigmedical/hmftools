package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.primitives.Doubles;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthEvidence;
import com.hartwig.hmftools.common.amber.BaseDepthFilter;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.amber.NormalHetrozygousFilter;
import com.hartwig.hmftools.common.amber.NormalHomozygousFilter;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAFEvidence;
import com.hartwig.hmftools.common.amber.TumorContamination;
import com.hartwig.hmftools.common.amber.TumorContaminationEvidence;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.region.GenomeRegion;

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

public class AmberApplication implements AutoCloseable {
    private static final Logger LOGGER = LogManager.getLogger(AmberApplication.class);

    private final AmberConfig config;
    private final ExecutorService executorService;
    private final Predicate<BaseDepth> homozygousFilter;
    private final Predicate<BaseDepth> heterozygousFilter;
    private final AmberPersistence persistence;

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

        final CommandLine cmd = createCommandLine(args, options);
        config = AmberConfig.createConfig(cmd);
        homozygousFilter = new NormalHomozygousFilter();
        heterozygousFilter = new NormalHetrozygousFilter(config.minHetAfPercent(), config.maxHetAfPercent());
        persistence = new AmberPersistence(config);

        final File outputDir = new File(config.outputDirectory());
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + config.outputDirectory());
        }

        if (!new File(config.bedFilePath()).exists()) {
            throw new IOException("Unable to locate bed file " + config.bedFilePath());
        }

        if (!new File(config.tumorBamPath()).exists()) {
            throw new IOException("Unable to locate tumor bam file " + config.tumorBamPath());
        }

        if (!new File(config.referenceBamPath()).exists()) {
            throw new IOException("Unable to locate reference bam file " + config.referenceBamPath());
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        executorService = Executors.newFixedThreadPool(config.threadCount(), namedThreadFactory);
    }

    private void run() throws InterruptedException, ExecutionException, IOException {
        final SamReaderFactory readerFactory = SamReaderFactory.make();

        final ListMultimap<Chromosome, BaseDepth> allNormal = normalDepth(readerFactory);
        final ListMultimap<Chromosome, BaseDepth> homNormal = Multimaps.filterEntries(allNormal, homozygousFilter);
        final ListMultimap<Chromosome, BaseDepth> hetNormal = Multimaps.filterEntries(allNormal, heterozygousFilter);
        final ListMultimap<Chromosome, TumorBAF> tumorBAFMap = tumorBAF(readerFactory, hetNormal);

        final List<TumorBAF> tumorBAFList =  tumorBAFMap.values().stream().sorted().collect(Collectors.toList());
        final List<AmberBAF> amberBAFList = tumorBAFList
                .stream()
                .map(AmberBAF::create)
                .filter(AmberApplication::isValid)
                .collect(Collectors.toList());

        final ListMultimap<Chromosome, TumorContamination> tumorContamination = contamination(readerFactory, homNormal);
        final List<TumorContamination> contaminationList = Lists.newArrayList(tumorContamination.values());

        persistence.persisQC(amberBAFList, contaminationList);
        persistence.persistVersionInfo();
        persistence.persistContamination(contaminationList);
        persistence.persistTumorBAF(tumorBAFList);
        persistence.persistAmberBAF(amberBAFList);
    }

    @NotNull
    private ListMultimap<Chromosome, BaseDepth> normalDepth(final SamReaderFactory readerFactory)
            throws IOException, InterruptedException, ExecutionException {
        LOGGER.info("Loading bed file {}", config.bedFilePath());
        final SortedSetMultimap<String, GenomeRegion> bedRegionsSortedSet = BEDFileLoader.fromBedFile(config.bedFilePath());
        final int partitionSize = Math.max(config.minPartition(), bedRegionsSortedSet.size() / config.threadCount());

        LOGGER.info("Processing {} potential sites in reference bam {}", bedRegionsSortedSet.values().size(), config.referenceBamPath());
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<BaseDepthEvidence>> futures = Lists.newArrayList();
        for (final String contig : bedRegionsSortedSet.keySet()) {
            for (final List<GenomeRegion> inner : Lists.partition(Lists.newArrayList(bedRegionsSortedSet.get(contig)), partitionSize)) {
                final BaseDepthEvidence evidence = new BaseDepthEvidence(config.typicalReadDepth(),
                        config.minMappingQuality(),
                        config.minBaseQuality(),
                        contig,
                        config.referenceBamPath(),
                        readerFactory,
                        inner);
                futures.add(executorService.submit(completion.task(evidence)));
            }
        }

        final ListMultimap<Chromosome, ModifiableBaseDepth> normalEvidence = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        final ListMultimap<Chromosome, BaseDepth> normalBafs = ArrayListMultimap.create();
        final Predicate<BaseDepth> depthFilter = new BaseDepthFilter(config.minDepthPercent(), config.maxDepthPercent(), normalEvidence);
        try (final RefEnricher refEnricher = new RefEnricher(config.refGenomePath())) {
            for (final Chromosome chromosome : normalEvidence.keySet()) {
                final List<BaseDepth> normalHetLocations = normalEvidence.get(chromosome)
                        .stream()
                        .filter(x -> x.indelCount() == 0)
                        .filter(depthFilter)
                        .map(refEnricher::enrich)
                        .filter(x -> heterozygousFilter.test(x) || homozygousFilter.test(x))
                        .collect(Collectors.toList());

                normalBafs.putAll(chromosome, normalHetLocations);
            }
        }

        return normalBafs;
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
        return Doubles.isFinite(baf.tumorBAF()) & Doubles.isFinite(baf.normalBAF());
    }

    @NotNull
    private static <T> List<T> getFuture(@NotNull final List<Future<T>> futures) throws ExecutionException, InterruptedException {
        final List<T> result = Lists.newArrayList();
        for (Future<T> chromosomeBAFEvidenceFuture : futures) {
            result.add(chromosomeBAFEvidenceFuture.get());
        }
        return result;
    }

    @Override
    public void close() {
        executorService.shutdown();
        LOGGER.info("Complete");
    }
}
