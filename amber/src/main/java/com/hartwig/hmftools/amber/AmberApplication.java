package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
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
import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.primitives.Doubles;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberVCF;
import com.hartwig.hmftools.common.amber.ImmutableAmberBAF;
import com.hartwig.hmftools.common.amber.ModifiableNormalBAF;
import com.hartwig.hmftools.common.amber.NormalBAF;
import com.hartwig.hmftools.common.amber.NormalBAFEvidence;
import com.hartwig.hmftools.common.amber.NormalDepthFilter;
import com.hartwig.hmftools.common.amber.NormalHetrozygousFilter;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAFEvidence;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFactory;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.version.VersionInfo;

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

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException {
        final Options options = AmberConfig.createOptions();
        try (final AmberApplication application = new AmberApplication(options, args)) {
            final List<TumorBAF> bafs = application.createBAFs();
            application.persist(bafs);
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

        final File outputDir = new File(config.outputDirectory());
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + config.outputDirectory());
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        executorService = Executors.newFixedThreadPool(config.threadCount(), namedThreadFactory);
    }

    @NotNull
    private List<TumorBAF> createBAFs() throws InterruptedException, ExecutionException, IOException {
        final SamReaderFactory readerFactory = SamReaderFactory.make();
        final ListMultimap<Chromosome, NormalBAF> normalBAFMap = normal(readerFactory);
        final Multimap<Chromosome, TumorBAF> tumorBAFMap = tumor(readerFactory, normalBAFMap);

        final List<TumorBAF> result = Lists.newArrayList(tumorBAFMap.values());
        Collections.sort(result);

        return result;
    }

    private void persist(@NotNull final List<TumorBAF> tumorBAFList) throws IOException, InterruptedException {
        LOGGER.info("Writing output to {}", config.outputDirectory());
        final String outputVcf = config.outputDirectory() + File.separator + config.tumor() + ".amber.vcf.gz";
        new AmberVCF(config.normal(), config.tumor()).write(outputVcf, tumorBAFList);

        final List<AmberBAF> result =
                tumorBAFList.stream().map(AmberApplication::create).filter(AmberApplication::isValid).collect(Collectors.toList());

        final AmberQC qcStats = AmberQCFactory.create(result);
        final String qcFilename = AmberQCFile.generateFilename(config.outputDirectory(), config.tumor());

        final String filename = AmberBAFFile.generateAmberFilename(config.outputDirectory(), config.tumor());
        AmberBAFFile.write(filename, result);
        AmberQCFile.write(qcFilename, qcStats);

        final VersionInfo versionInfo = new VersionInfo("amber.version");
        versionInfo.write(config.outputDirectory());

        LOGGER.info("Applying pcf segmentation");
        new BAFSegmentation(config.outputDirectory()).applySegmentation(config.tumor());
    }

    @NotNull
    private ListMultimap<Chromosome, TumorBAF> tumor(@NotNull final SamReaderFactory readerFactory,
            @NotNull final ListMultimap<Chromosome, NormalBAF> normalBafs) throws ExecutionException, InterruptedException {
        final int partitionSize = Math.max(config.minPartition(), normalBafs.values().size() / config.threadCount());

        LOGGER.info("Processing tumor bam {}", config.tumorBamPath());
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<TumorBAFEvidence>> futures = Lists.newArrayList();
        for (final Chromosome chromosome : normalBafs.keySet()) {
            for (final List<NormalBAF> chromosomeBafPoints : Lists.partition(normalBafs.get(chromosome), partitionSize)) {
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
    private ListMultimap<Chromosome, NormalBAF> normal(final SamReaderFactory readerFactory)
            throws IOException, InterruptedException, ExecutionException {
        LOGGER.info("Loading bed file {}", config.bedFilePath());
        final SortedSetMultimap<String, GenomeRegion> bedRegionsSortedSet = BEDFileLoader.fromBedFile(config.bedFilePath());
        final int partitionSize = Math.max(config.minPartition(), bedRegionsSortedSet.size() / config.threadCount());

        LOGGER.info("Processing reference bam {}", config.referenceBamPath(), partitionSize);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<NormalBAFEvidence>> futures = Lists.newArrayList();
        for (final String contig : bedRegionsSortedSet.keySet()) {
            for (final List<GenomeRegion> inner : Lists.partition(Lists.newArrayList(bedRegionsSortedSet.get(contig)), partitionSize)) {
                final NormalBAFEvidence evidence = new NormalBAFEvidence(config.typicalReadDepth(),
                        config.minMappingQuality(),
                        config.minBaseQuality(),
                        contig,
                        config.referenceBamPath(),
                        readerFactory,
                        inner);
                futures.add(executorService.submit(completion.task(evidence)));
            }
        }

        final ListMultimap<Chromosome, ModifiableNormalBAF> normalEvidence = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        final ListMultimap<Chromosome, NormalBAF> normalBafs = ArrayListMultimap.create();
        final Predicate<NormalBAF> depthFilter = new NormalDepthFilter(config.minDepthPercent(), config.maxDepthPercent(), normalEvidence);
        final RefEnricher refEnricher = new RefEnricher(config.refGenomePath());
        final Predicate<NormalBAF> hetFilter = new NormalHetrozygousFilter(config.minHetAfPercent(), config.maxHetAfPercent());
        for (final Chromosome chromosome : normalEvidence.keySet()) {
            final List<NormalBAF> normalHetLocations = normalEvidence.get(chromosome)
                    .stream()
                    .filter(x -> x.indelCount() == 0)
                    .filter(depthFilter)
                    .map(refEnricher::enrich)
                    .filter(hetFilter)
                    .collect(Collectors.toList());

            normalBafs.putAll(chromosome, normalHetLocations);
        }

        LOGGER.info("Identified {} heterozygous sites", normalBafs.values().size());
        return normalBafs;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static AmberBAF create(@NotNull final TumorBAF tumor) {
        int tumorAltCount = tumor.tumorAltSupport();
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumor.tumorRefSupport());
        int normalAltCount = tumor.normalAltSupport();
        double normalBaf = normalAltCount / (double) (normalAltCount + tumor.normalRefSupport());
        return ImmutableAmberBAF.builder()
                .from(tumor)
                .normalDepth(tumor.normalReadDepth())
                .tumorDepth(tumor.tumorReadDepth())
                .normalBAF(normalBaf)
                .tumorBAF(tumorBaf)
                .build();
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
