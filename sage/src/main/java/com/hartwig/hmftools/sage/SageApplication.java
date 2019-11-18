package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.pipeline.SlicePipeline;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class SageApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final SageVCF vcf;
    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final ListMultimap<Chromosome, VariantHotspot> hotspots;
    private final ListMultimap<Chromosome, GenomeRegion> panel = ArrayListMultimap.create();
    private final SamSlicerFactory samSlicerFactory;

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException {
        final Options options = SageConfig.createOptions();
        try (final SageApplication application = new SageApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageApplication", options);
            System.exit(1);
        }
    }

    private SageApplication(final Options options, final String... args) throws IOException, ParseException {
        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(cmd);

        if (!config.panel().isEmpty()) {
            LOGGER.info("Reading gene panel bed file: {}", config.panel());
            SortedSetMultimap<String, GenomeRegion> bed = BEDFileLoader.fromBedFile(config.panel());
            for (String contig : bed.keySet()) {
                if (HumanChromosome.contains(contig)) {
                    panel.putAll(HumanChromosome.fromString(contig), bed.get(contig));
                }
            }
        }

        if (!config.hotspots().isEmpty()) {
            LOGGER.info("Reading hotspot vcf: {}", config.hotspots());
            hotspots = VariantHotspotFile.readFromVCF(config.hotspots());
        } else {
            hotspots = ArrayListMultimap.create();
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        vcf = new SageVCF(refGenome, config);
        samSlicerFactory = new SamSlicerFactory(config, panel);
        LOGGER.info("Writing to file {}", config.outputFile());
    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        long timeStamp = System.currentTimeMillis();
        final List<Future<ChromosomePipeline>> contigContexts = Lists.newArrayList();

        SAMSequenceDictionary dictionary = dictionary();
        for (final SAMSequenceRecord samSequenceRecord : dictionary.getSequences()) {
            final String contig = samSequenceRecord.getSequenceName();
            if (HumanChromosome.contains(contig)) {
                int maxPosition = samSequenceRecord.getSequenceLength();
                contigContexts.add(runChromosome(contig, config.regionSliceSize(), maxPosition));
            }
        }

        //        contigContexts.add(runChromosome("17", config.regionSliceSize(), 4_000_000));
        //        contigContexts.add(runChromosome("17", config.regionSliceSize(), dictionary().getSequence("17").getSequenceLength()));
        //                contigContexts.add(runSingleRegion("17", 6133723, 6133723));
        //                        contigContexts.add(runSingleRegion("1", 696644, 696644));
        //                contigContexts.add(runSingleRegion("17", 3028422, 3028422));
        //        contigContexts.add(runSingleRegion("17", 19_465_877, 19465877));
        //                contigContexts.add(runSingleRegion("17", 1743210, 1743211));
        //        contigContexts.add(runSingleRegion("17", 22_260_001, 23_262_000));
        //        contigContexts.add(runSingleRegion("17", 25_282_540, 34000000));
        //                contigContexts.add(runSingleRegion("1", 159946533, 159946943));
        //        contigContexts.add(runSingleRegion("17", 37_000_000, 38_000_000));
        //        contigContexts.add(runSingleRegion("17", 42_796_634, 42796634));
        //        contigContexts.add(runSingleRegion("17", 47_414_327, 47414327));
        //        contigContexts.add(runSingleRegion("17", 55_639_513, 55639513));
        //        contigContexts.add(runSingleRegion("17", 72_558_371, 72558371));
        //                        contigContexts.add(runSingleRegion("17", 22163006, 25363006));
        //                contigContexts.add(runSingleRegion("4", 943940, 943950));
        //                contigContexts.add(runSingleRegion("5", 68706692, 68706892));
        //                contigContexts.add(runSingleRegion("1", 51486496, 51486496));
        //                        contigContexts.add(runSingleRegion("22", 29453442, 29453442));
        //                                contigContexts.add(runSingleRegion("17", 22_200_000, 22_300_000));
        //                                contigContexts.add(runSingleRegion("10", 42350000, 42450000));
        //                                contigContexts.add(runSingleRegion("4", 49000000, 50000000));

        for (Future<ChromosomePipeline> contigContext : contigContexts) {
            final ChromosomePipeline chromosomePipeline = contigContext.get();
            LOGGER.info("Writing chromosome {} ", chromosomePipeline.chromosome());
            for (VariantContext context : chromosomePipeline.variantContexts()) {
                vcf.write(context);
            }
        }

        long timeTaken = System.currentTimeMillis() - timeStamp;
        LOGGER.info("Completed in {} seconds", timeTaken / 1000);
    }

    private SAMSequenceDictionary dictionary() throws IOException {
        SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(config.referenceBam()));
        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

    @NotNull
    private Future<ChromosomePipeline> runChromosome(@NotNull final String contig, int regionSliceSize, int maxPosition) {
        final Chromosome chromosome = HumanChromosome.fromString(contig);
        final SageVariantFactory variantFactory =
                new SageVariantFactory(chromosome, config.filter(), hotspots.get(chromosome), panel.get(chromosome));
        final ChromosomePipeline chromosomePipeline = new ChromosomePipeline(contig, config, refGenome, variantFactory);
        for (int i = 0; ; i++) {
            int start = 1 + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize, maxPosition);
            chromosomePipeline.accept(runRegion(contig, start, end));

            if (end >= maxPosition) {
                break;
            }
        }

        return chromosomePipeline.submit();
    }

    @NotNull
    private Future<ChromosomePipeline> runSingleRegion(@NotNull final String contig, int start, int end) {
        final Chromosome chromosome = HumanChromosome.fromString(contig);
        final SageVariantFactory variantFactory =
                new SageVariantFactory(chromosome, config.filter(), hotspots.get(chromosome), panel.get(chromosome));
        final ChromosomePipeline chromosomePipeline = new ChromosomePipeline(contig, config, refGenome, variantFactory);

        chromosomePipeline.accept(runRegion(contig, start, end));
        return chromosomePipeline.submit();

    }

    @NotNull
    private CompletableFuture<List<SageVariant>> runRegion(@NotNull final String contig, int start, int end) {
        final Chromosome chromosome = HumanChromosome.fromString(contig);
        final SageVariantFactory variantFactory =
                new SageVariantFactory(chromosome, config.filter(), hotspots.get(chromosome), panel.get(chromosome));

        final GenomeRegion region = GenomeRegions.create(contig, start, end);
        final SlicePipeline regionPipeline =
                new SlicePipeline(region, config, executorService, refGenome, variantFactory, samSlicerFactory);
        return regionPipeline.submit();

    }

    @Override
    public void close() throws IOException {
        vcf.close();
        refGenome.close();
        executorService.shutdown();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
