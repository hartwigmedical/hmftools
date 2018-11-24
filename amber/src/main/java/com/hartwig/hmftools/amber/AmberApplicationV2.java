package com.hartwig.hmftools.amber;

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

    static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;
    static final double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    static final double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;

    private static final String THREADS = "threads";
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String BED_FILE = "bed";
    private static final String REF_GENOME = "ref_genome";

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

        final File outputDir = new File(outputDirectory);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + outputDirectory);
        }

        final int threadCount = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : 1;
        final SamReaderFactory readerFactory = SamReaderFactory.make();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        final ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);

        final ListMultimap<Chromosome, NormalBAF> normalBAFMap =
                normal(readerFactory, executorService, bedFilePath, refGenomePath, referenceBamPath);
        final Multimap<Chromosome, TumorBAF> tumorBAFMap = tumor(tumorBamPath, readerFactory, executorService, normalBAFMap);

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
            @NotNull final ListMultimap<Chromosome, NormalBAF> normalBafs) throws ExecutionException, InterruptedException {
        LOGGER.info("Processing tumor bam file {}", tumorBamPath);
        final List<Future<TumorEvidence>> futures = Lists.newArrayList();
        for (final Chromosome chromosome : normalBafs.keySet()) {
            final List<NormalBAF> chromosomeBafPoints = normalBafs.get(chromosome);
            if (!chromosomeBafPoints.isEmpty()) {
                final String contig = chromosomeBafPoints.get(0).chromosome();
                futures.add(executorService.submit(new TumorEvidence(contig, tumorBamPath, readerFactory, chromosomeBafPoints)));
            }
        }

        final ListMultimap<Chromosome, TumorBAF> result = ArrayListMultimap.create();
        for (Future<TumorEvidence> future : futures) {
            final TumorEvidence evidence = future.get();
            result.putAll(HumanChromosome.fromString(evidence.contig()), evidence.evidence());
        }

        return result;
    }

    @NotNull
    private static ListMultimap<Chromosome, NormalBAF> normal(final SamReaderFactory readerFactory, final ExecutorService executorService,
            final String bedPath, final String refGenomePath, final String referenceBamPath)
            throws IOException, InterruptedException, ExecutionException {
        LOGGER.info("Loading bed file {}", bedPath);
        final SortedSetMultimap<String, GenomeRegion> bedRegionsSortedSet = BEDFileLoader.fromBedFile(bedPath);

        LOGGER.info("Processing reference bam file {}", referenceBamPath);
        final List<Future<NormalEvidence>> futures = Lists.newArrayList();
        for (final String contig : bedRegionsSortedSet.keySet()) {
            final NormalEvidence chromosome =
                    new NormalEvidence(contig, referenceBamPath, readerFactory, new ArrayList<>(bedRegionsSortedSet.get(contig)));
            futures.add(executorService.submit(chromosome));
        }
        final ListMultimap<Chromosome, ModifiableNormalBAF> normalEvidence = ArrayListMultimap.create();
        for (Future<NormalEvidence> chromosomeBAFEvidenceFuture : futures) {
            NormalEvidence evidence = chromosomeBAFEvidenceFuture.get();
            normalEvidence.putAll(HumanChromosome.fromString(evidence.contig()), evidence.getEvidence());
        }

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

        LOGGER.info("Identified {} heterozygous sites in reference bam", normalBafs.values().size());
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
        options.addOption(THREADS, true, "Number of threads. Default 1.");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample.");
        options.addOption(TUMOR_BAM, true, "Tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(BED_FILE, true, "Baf locations bed file.");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
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
}
