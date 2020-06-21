package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.function.BiFunction;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.quality.QualityRecalibration;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationFile;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.vcf.SageVCF;

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
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final SageVCF vcf;
    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;

    private final ListMultimap<Chromosome, GenomeRegion> panel;
    private final ListMultimap<Chromosome, VariantHotspot> hotspots;
    private final ListMultimap<Chromosome, GenomeRegion> highConfidence;

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
        final VersionInfo version = new VersionInfo("sage.version");
        LOGGER.info("SAGE version: {}", version.version());

        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(version.version(), cmd);

        hotspots = readHotspots();
        panel = panelWithHotspots(hotspots);
        highConfidence = readPanel(config.highConfidenceBed());

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));

        vcf = new SageVCF(refGenome, config);
        LOGGER.info("Writing to file: {}", config.outputFile());
    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        long timeStamp = System.currentTimeMillis();

        final Map<String, QualityRecalibrationMap> recalibrationMap = qualityRecalibration();
        final SAMSequenceDictionary dictionary = dictionary();
        for (final SAMSequenceRecord samSequenceRecord : dictionary.getSequences()) {
            final String contig = samSequenceRecord.getSequenceName();
            if (config.chromosomes().isEmpty() || config.chromosomes().contains(contig)) {
                if (HumanChromosome.contains(contig) || MitochondrialChromosome.contains(contig)) {
                    try (final ChromosomePipeline pipeline = createChromosomePipeline(contig, recalibrationMap)) {
                        pipeline.process();
                    }
                    System.gc();
                }
            }
        }

//                createChromosomePipeline("10", recalibrationMap).process(130404941, 130405950);

        long timeTaken = System.currentTimeMillis() - timeStamp;
        LOGGER.info("Completed in {} seconds", timeTaken / 1000);
    }

    private SAMSequenceDictionary dictionary() throws IOException {
        final String bam = config.referenceBam().isEmpty() ? config.tumorBam().get(0) : config.referenceBam().get(0);
        SamReader tumorReader = SamReaderFactory.makeDefault().referenceSource(new ReferenceSource(refGenome)).open(new File(bam));
        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

    private ChromosomePipeline createChromosomePipeline(@NotNull final String contig,
            Map<String, QualityRecalibrationMap> qualityRecalibrationMap) throws IOException {
        final Chromosome chromosome =
                HumanChromosome.contains(contig) ? HumanChromosome.fromString(contig) : MitochondrialChromosome.fromString(contig);
        return new ChromosomePipeline(contig,
                config,
                executorService,
                hotspots.get(chromosome),
                panel.get(chromosome),
                highConfidence.get(chromosome),
                qualityRecalibrationMap,
                vcf::write);
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

    @NotNull
    private ListMultimap<Chromosome, VariantHotspot> readHotspots() throws IOException {
        if (!config.hotspots().isEmpty()) {
            LOGGER.info("Reading hotspot vcf: {}", config.hotspots());
            return VariantHotspotFile.readFromVCF(config.hotspots());
        } else {
            return ArrayListMultimap.create();
        }
    }

    @NotNull
    private ListMultimap<Chromosome, GenomeRegion> panelWithHotspots(@NotNull final ListMultimap<Chromosome, VariantHotspot> hotspots)
            throws IOException {
        final ListMultimap<Chromosome, GenomeRegion> initialPanel = readPanel(config.panelBed());
        final ListMultimap<Chromosome, GenomeRegion> result = ArrayListMultimap.create();

        for (HumanChromosome chromosome : HumanChromosome.values()) {
            final GenomeRegions builder = new GenomeRegions(chromosome.toString());
            if (initialPanel.containsKey(chromosome)) {
                initialPanel.get(chromosome).forEach(x -> builder.addRegion(x.start(), x.end()));
            }
            if (hotspots.containsKey(chromosome)) {
                hotspots.get(chromosome).forEach(x -> builder.addPosition(x.position()));
            }

            result.putAll(chromosome, builder.build());
        }

        return result;
    }

    @NotNull
    private static ListMultimap<Chromosome, GenomeRegion> readPanel(@NotNull final String panelBed) throws IOException {
        final ListMultimap<Chromosome, GenomeRegion> panel = ArrayListMultimap.create();
        if (!panelBed.isEmpty()) {
            LOGGER.info("Reading bed file: {}", panelBed);
            SortedSetMultimap<String, GenomeRegion> bed = BEDFileLoader.fromBedFile(panelBed);
            for (String contig : bed.keySet()) {
                if (HumanChromosome.contains(contig)) {
                    panel.putAll(HumanChromosome.fromString(contig), bed.get(contig));
                }
            }
        }

        return panel;
    }

    @NotNull
    private Map<String, QualityRecalibrationMap> qualityRecalibration() {
        final BaseQualityRecalibrationConfig bqrConfig = config.baseQualityRecalibrationConfig();

        if (!bqrConfig.enabled()) {
            return disableQualityRecalibration();
        }

        final Map<String, QualityRecalibrationMap> result = Maps.newHashMap();
        LOGGER.info("Beginning quality recalibration");

        final QualityRecalibration qualityRecalibration = new QualityRecalibration(bqrConfig, executorService, refGenome);
        final List<CompletableFuture<Void>> done = Lists.newArrayList();

        final BiFunction<String, String, CompletableFuture<Void>> processSample =
                (sample, sampleBam) -> qualityRecalibration.qualityRecalibrationRecords(sampleBam).thenAccept(records -> {
                    try {

                        final String tsvFile = config.baseQualityRecalibrationFile(sample);
                        QualityRecalibrationFile.write(tsvFile, records);
                        result.put(sample, new QualityRecalibrationMap(records));
                        LOGGER.info("Writing base quality recalibration file: {}", tsvFile);
                        if (bqrConfig.plot()) {
                            RExecutor.executeFromClasspath("r/baseQualityRecalibrationPlot.R", tsvFile);
                        }
                    } catch (Exception e) {
                        throw new CompletionException(e);
                    }
                });

        for (int i = 0; i < config.reference().size(); i++) {
            done.add(processSample.apply(config.reference().get(i), config.referenceBam().get(i)));
        }

        for (int i = 0; i < config.tumor().size(); i++) {
            done.add(processSample.apply(config.tumor().get(i), config.tumorBam().get(i)));
        }

        // Wait for all tasks to be finished
        done.forEach(CompletableFuture::join);

        return result;
    }

    private Map<String, QualityRecalibrationMap> disableQualityRecalibration() {
        final Map<String, QualityRecalibrationMap> result = Maps.newHashMap();

        for (String sample : config.reference()) {
            result.put(sample, new QualityRecalibrationMap(Collections.emptyList()));
        }
        for (String sample : config.tumor()) {
            result.put(sample, new QualityRecalibrationMap(Collections.emptyList()));

        }
        return result;
    }

}