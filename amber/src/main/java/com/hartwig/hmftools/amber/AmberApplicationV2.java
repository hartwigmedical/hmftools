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

import com.google.common.base.Strings;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.primitives.Doubles;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.ImmutableAmberBAF;
import com.hartwig.hmftools.common.amber.ModifiableNormalBAF;
import com.hartwig.hmftools.common.amber.NormalBAF;
import com.hartwig.hmftools.common.amber.TumorBAF;
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

public class AmberApplicationV2 {
    private static final Logger LOGGER = LogManager.getLogger(AmberApplicationV2.class);

    private static final int MIN_PARITION = 10000;
    private static final int DEFAULT_THREADS = 1;
    private static final int DEFAULT_MIN_BASE_QUALITY = 13;
    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;
    private static final double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    private static final double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;
    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    private static final String TUMOR = "tumor";
    private static final String BED_FILE = "bed";
    private static final String THREADS = "threads";
    private static final String REFERENCE = "reference";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String REF_GENOME = "ref_genome";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String MIN_BASE_QUALITY = "min_base_quality";
    private static final String MIN_MAPPING_QUALITY = "min_mapping_quality";

    private static final String MIN_DEPTH_PERCENTAGE = "min_depth_percent";
    private static final String MAX_DEPTH_PERCENTAGE = "max_depth_percent";
    private static final String MIN_HET_AF_PERCENTAGE = "min_het_af_percent";
    private static final String MAX_HET_AF_PERCENTAGE = "max_het_af_percent";

    public static void main(final String... args) throws ParseException, IOException, InterruptedException, ExecutionException {
        final VersionInfo versionInfo = new VersionInfo("amber.version");
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String bedFilePath = cmd.getOptionValue(BED_FILE);
        final String tumorBamPath = cmd.getOptionValue(TUMOR_BAM);
        final String referenceBamPath = cmd.getOptionValue(REFERENCE_BAM);
        final String refGenomePath = cmd.getOptionValue(REF_GENOME);
        final String outputDirectory = cmd.getOptionValue(OUTPUT_DIR);
        final String normal = cmd.getOptionValue(REFERENCE);
        final String tumor = cmd.getOptionValue(TUMOR);

        if (bedFilePath == null || referenceBamPath == null || tumorBamPath == null || outputDirectory == null || normal == null
                || tumor == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("AbamberApplication", options);
            System.exit(1);
        }

        int minBaseQuality =
                cmd.hasOption(MIN_BASE_QUALITY) ? Integer.valueOf(cmd.getOptionValue(MIN_BASE_QUALITY)) : DEFAULT_MIN_BASE_QUALITY;

        final File outputDir = new File(outputDirectory);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + outputDirectory);
        }

        final int threadCount = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : DEFAULT_THREADS;
        final SamReaderFactory readerFactory = SamReaderFactory.make();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        final ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);

        final ListMultimap<Chromosome, NormalBAF> normalBAFMap =
                normal(readerFactory, executorService, bedFilePath, refGenomePath, referenceBamPath, minBaseQuality, threadCount);
        final Multimap<Chromosome, TumorBAF> tumorBAFMap =
                tumor(tumorBamPath, readerFactory, executorService, normalBAFMap, minBaseQuality, threadCount);

        final List<TumorBAF> tumorBAFList = Lists.newArrayList(tumorBAFMap.values());
        Collections.sort(tumorBAFList);

        final String outputVcf = outputDirectory + File.separator + tumor + ".amber.vcf.gz";
        LOGGER.info("Writing vcf output to {}", outputVcf);
        new AmberVCF(normal, tumor).write(outputVcf, tumorBAFMap.values());

        final List<AmberBAF> result =
                tumorBAFList.stream().map(AmberApplicationV2::create).filter(AmberApplicationV2::isValid).collect(Collectors.toList());

        LOGGER.info("Generating QC Stats");
        final AmberQC qcStats = AmberQCFactory.create(result);
        final String qcFilename = AmberQCFile.generateFilename(outputDirectory, tumor);

        final String filename = AmberBAFFile.generateAmberFilename(outputDirectory, tumor);
        LOGGER.info("Persisting file {}", filename);
        AmberBAFFile.write(filename, result);
        AmberQCFile.write(qcFilename, qcStats);
        versionInfo.write(outputDir.toString());

        LOGGER.info("Applying pcf segmentation");
        new BAFSegmentation(outputDirectory).applySegmentation(tumor);

        executorService.shutdown();
        LOGGER.info("Complete");
    }

    @NotNull
    private static ListMultimap<Chromosome, TumorBAF> tumor(@NotNull final String tumorBamPath,
            @NotNull final SamReaderFactory readerFactory, @NotNull final ExecutorService executorService,
            @NotNull final ListMultimap<Chromosome, NormalBAF> normalBafs, int minBaseQuality, int threads)
            throws ExecutionException, InterruptedException {
        final int partitionSize = Math.max(MIN_PARITION, normalBafs.values().size() / threads);

        LOGGER.info("Processing tumor bam {}", tumorBamPath);
        final List<Future<TumorEvidence>> futures = Lists.newArrayList();
        for (final Chromosome chromosome : normalBafs.keySet()) {
            for (final List<NormalBAF> chromosomeBafPoints : Lists.partition(normalBafs.get(chromosome), partitionSize)) {
                if (!chromosomeBafPoints.isEmpty()) {
                    final String contig = chromosomeBafPoints.get(0).chromosome();
                    futures.add(executorService.submit(new TumorEvidence(minBaseQuality,
                            contig,
                            tumorBamPath,
                            readerFactory,
                            chromosomeBafPoints)));
                }
            }
        }

        final ListMultimap<Chromosome, TumorBAF> result = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> result.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        return result;
    }

    @NotNull
    private static ListMultimap<Chromosome, NormalBAF> normal(final SamReaderFactory readerFactory, final ExecutorService executorService,
            final String bedPath, final String refGenomePath, final String referenceBamPath, int minBaseQuality, int threads)
            throws IOException, InterruptedException, ExecutionException {
        LOGGER.info("Loading bed file {}", bedPath);
        final SortedSetMultimap<String, GenomeRegion> bedRegionsSortedSet = BEDFileLoader.fromBedFile(bedPath);
        final int partitionSize = Math.max(MIN_PARITION, bedRegionsSortedSet.size() / threads);

        LOGGER.info("Processing reference bam {}", referenceBamPath, partitionSize);
        final List<Future<NormalEvidence>> futures = Lists.newArrayList();
        for (final String contig : bedRegionsSortedSet.keySet()) {
            for (final List<GenomeRegion> inner : Lists.partition(Lists.newArrayList(bedRegionsSortedSet.get(contig)), partitionSize)) {
                final NormalEvidence chromosome = new NormalEvidence(minBaseQuality, contig, referenceBamPath, readerFactory, inner);
                futures.add(executorService.submit(chromosome));
            }
        }

        final ListMultimap<Chromosome, ModifiableNormalBAF> normalEvidence = ArrayListMultimap.create();
        getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));

        final ListMultimap<Chromosome, NormalBAF> normalBafs = ArrayListMultimap.create();
        final Predicate<NormalBAF> depthFilter =
                new DepthFilter(DEFAULT_MIN_DEPTH_PERCENTAGE, DEFAULT_MAX_DEPTH_PERCENTAGE, normalEvidence);
        final RefEnricher refEnricher = new RefEnricher(refGenomePath);
        final Predicate<NormalBAF> hetFilter = new HetrozygousFilter(DEFAULT_MIN_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);
        for (final Chromosome chromosome : normalEvidence.keySet()) {
            final List<NormalBAF> normalHetLocations = normalEvidence.get(chromosome)
                    .stream()
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
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample.");
        options.addOption(TUMOR_BAM, true, "Tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(BED_FILE, true, "Baf locations bed file.");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(MIN_BASE_QUALITY, true, "Minimum base quality for a base to be considered [" + DEFAULT_MIN_BASE_QUALITY + "]");
        options.addOption(MIN_MAPPING_QUALITY,
                true,
                "Minimum mapping quality for an alignment to be used [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(MIN_HET_AF_PERCENTAGE, true, "Min heterozygous AF% [" + DEFAULT_MIN_HET_AF_PERCENTAGE + "]");
        options.addOption(MAX_HET_AF_PERCENTAGE, true, "Max heterozygous AF% [" + DEFAULT_MAX_HET_AF_PERCENTAGE + "]");
        options.addOption(MIN_DEPTH_PERCENTAGE, true, "Max percentage of median depth [" + DEFAULT_MIN_DEPTH_PERCENTAGE + "]");
        options.addOption(MAX_DEPTH_PERCENTAGE, true, "Min percentage of median depth [" + DEFAULT_MAX_DEPTH_PERCENTAGE + "]");
        return options;
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
        int complete = 0;
        double previousPercentComplete = 0;
        for (Future<T> chromosomeBAFEvidenceFuture : futures) {
            T evidence = chromosomeBAFEvidenceFuture.get();
            double percentComplete = ((double) ++complete) / futures.size();
            if (percentComplete > previousPercentComplete + 0.1) {
                LOGGER.info("{}", complete(percentComplete));
                previousPercentComplete = percentComplete;
            }
            result.add(evidence);
        }
        return result;
    }

    @NotNull
    private static String complete(double percent) {
        int roundedPercent = (int) Math.round(percent * 100);
        int hashCount = Math.min(20, roundedPercent / 5);
        int gapCount = Math.max(0, 20 - hashCount);

        return "  [" + Strings.repeat("#", hashCount) + Strings.repeat(" ", gapCount) + "] " + roundedPercent + "% complete";
    }
}
