package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsValidation;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneDepthFile;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationSupplier;
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

public class SageApplication implements AutoCloseable
{

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final SageVCF vcf;
    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final QualityRecalibrationSupplier qualityRecalibrationSupplier;

    private final ListMultimap<Chromosome, NamedBed> coveragePanel;
    private final ListMultimap<Chromosome, GenomeRegion> panelWithHotspots;
    private final ListMultimap<Chromosome, VariantHotspot> hotspots;
    private final ListMultimap<Chromosome, GenomeRegion> highConfidence;

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SageConfig.createSageOptions();
        try(final SageApplication application = new SageApplication(options, args))
        {
            application.run();
        } catch(ParseException e)
        {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageApplication", options);
            System.exit(1);
        }
    }

    private SageApplication(final Options options, final String... args) throws IOException, ParseException
    {
        final VersionInfo version = new VersionInfo("sage.version");
        LOGGER.info("SAGE version: {}", version.version());

        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(false, version.version(), cmd);
        coveragePanel = readNamedBed(config.coverageBed());
        final ListMultimap<Chromosome, GenomeRegion> panelWithoutHotspots = readUnnamedBed(config.panelBed());
        hotspots = readHotspots();
        panelWithHotspots = panelWithHotspots(panelWithoutHotspots, hotspots);
        highConfidence = readUnnamedBed(config.highConfidenceBed());

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        qualityRecalibrationSupplier = new QualityRecalibrationSupplier(executorService, refGenome, config);

        vcf = new SageVCF(refGenome, config);
        LOGGER.info("Writing to file: {}", config.outputFile());

        // Validate Coverage Bed
        if(config.panelOnly() && !coveragePanel.isEmpty())
        {
            if(!GenomeRegionsValidation.isSubset(panelWithHotspots.values(), coveragePanel.values()))
            {
                throw new IOException("Coverage bed must be a subset of panel bed when running in panel only mode");
            }
        }
    }

    @NotNull
    private Coverage createCoverage()
    {
        Set<String> samples = Sets.newHashSet();
        if(!config.coverageBed().isEmpty())
        {
            samples.addAll(config.tumor());
        }
        return new Coverage(samples, coveragePanel.values());
    }

    private void run() throws InterruptedException, ExecutionException, IOException
    {
        long timeStamp = System.currentTimeMillis();
        final Coverage coverage = createCoverage();

        final Map<String, QualityRecalibrationMap> recalibrationMap = qualityRecalibrationSupplier.get();
        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String contig = samSequenceRecord.getSequenceName();
            if(config.chromosomes().isEmpty() || config.chromosomes().contains(contig))
            {
                if(HumanChromosome.contains(contig) || MitochondrialChromosome.contains(contig))
                {
                    try(final ChromosomePipeline pipeline = createChromosomePipeline(contig, coverage, recalibrationMap))
                    {
                        pipeline.process();
                    }
                    System.gc();
                }
            }
        }

        //                createChromosomePipeline("10", recalibrationMap).process(130404941, 130405950);

        // Write out coverage
        for(String sample : coverage.samples())
        {
            String filename = config.geneCoverageFile(sample);
            GeneDepthFile.write(filename, coverage.depth(sample));
        }

        long timeTaken = System.currentTimeMillis() - timeStamp;
        LOGGER.info("Completed in {} seconds", timeTaken / 1000);
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = config.referenceBam().isEmpty() ? config.tumorBam().get(0) : config.referenceBam().get(0);
        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(config.validationStringency())
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bam));
        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

    private ChromosomePipeline createChromosomePipeline(@NotNull final String contig, @NotNull final Coverage coverage,
            @NotNull Map<String, QualityRecalibrationMap> qualityRecalibrationMap) throws IOException
    {
        final Chromosome chromosome =
                HumanChromosome.contains(contig) ? HumanChromosome.fromString(contig) : MitochondrialChromosome.fromString(contig);
        return new ChromosomePipeline(contig,
                config,
                executorService,
                hotspots.get(chromosome),
                panelWithHotspots.get(chromosome),
                highConfidence.get(chromosome),
                qualityRecalibrationMap,
                coverage,
                vcf::write);
    }

    @Override
    public void close() throws IOException
    {
        vcf.close();
        refGenome.close();
        executorService.shutdown();
    }

    @NotNull
    static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private ListMultimap<Chromosome, VariantHotspot> readHotspots() throws IOException
    {
        if(!config.hotspots().isEmpty())
        {
            LOGGER.info("Reading hotspot vcf: {}", config.hotspots());
            return VariantHotspotFile.readFromVCF(config.hotspots());
        }
        else
        {
            return ArrayListMultimap.create();
        }
    }

    @NotNull
    private ListMultimap<Chromosome, GenomeRegion> panelWithHotspots(final ListMultimap<Chromosome, GenomeRegion> panelWithoutHotspots,
            @NotNull final ListMultimap<Chromosome, VariantHotspot> hotspots)
    {
        final ListMultimap<Chromosome, GenomeRegion> result = ArrayListMultimap.create();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            final GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
            if(panelWithoutHotspots.containsKey(chromosome))
            {
                panelWithoutHotspots.get(chromosome).forEach(builder::addRegion);
            }
            if(hotspots.containsKey(chromosome))
            {
                hotspots.get(chromosome).forEach(builder::addPosition);
            }

            result.putAll(chromosome, builder.build());
        }

        return result;
    }

    @NotNull
    private static ListMultimap<Chromosome, NamedBed> readNamedBed(@NotNull final String panelBed) throws IOException
    {
        final ListMultimap<Chromosome, NamedBed> panel = ArrayListMultimap.create();
        if(!panelBed.isEmpty())
        {
            LOGGER.info("Reading bed file: {}", panelBed);
            for(NamedBed bed : NamedBedFile.readBedFile(panelBed))
            {
                if(HumanChromosome.contains(bed.chromosome()))
                {
                    panel.put(HumanChromosome.fromString(bed.chromosome()), bed);
                }
            }
        }

        return panel;
    }

    @NotNull
    private static ListMultimap<Chromosome, GenomeRegion> readUnnamedBed(@NotNull final String panelBed) throws IOException
    {
        final ListMultimap<Chromosome, GenomeRegion> panel = ArrayListMultimap.create();
        if(!panelBed.isEmpty())
        {
            LOGGER.info("Reading bed file: {}", panelBed);
            SortedSetMultimap<String, GenomeRegion> bed = BEDFileLoader.fromBedFile(panelBed);
            for(String contig : bed.keySet())
            {
                if(HumanChromosome.contains(contig))
                {
                    panel.putAll(HumanChromosome.fromString(contig), bed.get(contig));
                }
            }
        }

        return panel;
    }

}