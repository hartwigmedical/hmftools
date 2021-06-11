package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

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
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageApplication implements AutoCloseable
{
    private final SageVCF mVcfFile;
    private final SageConfig mConfig;
    private final ExecutorService mExecutorService;
    private final IndexedFastaSequenceFile mRefGenome;
    private final QualityRecalibrationSupplier mQualityRecalibrationSupplier;

    private final ListMultimap<Chromosome, NamedBed> mCoveragePanel;
    private final ListMultimap<Chromosome, GenomeRegion> mPanelWithHotspots;
    private final ListMultimap<Chromosome, VariantHotspot> mHotspots;
    private final ListMultimap<Chromosome, GenomeRegion> mHighConfidence;

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SageConfig.createSageOptions();

        try(final SageApplication application = new SageApplication(options, args))
        {
            application.run();
        }
        catch(ParseException e)
        {
            SG_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageApplication", options);
            System.exit(1);
        }
    }

    private SageApplication(final Options options, final String... args) throws IOException, ParseException
    {
        final VersionInfo version = new VersionInfo("sage.version");
        SG_LOGGER.info("SAGE version: {}", version.version());

        final CommandLine cmd = createCommandLine(args, options);
        this.mConfig = SageConfig.createConfig(false, version.version(), cmd);
        mCoveragePanel = readNamedBed(mConfig.coverageBed());
        final ListMultimap<Chromosome, GenomeRegion> panelWithoutHotspots = readUnnamedBed(mConfig.panelBed());
        mHotspots = readHotspots();
        mPanelWithHotspots = panelWithHotspots(panelWithoutHotspots, mHotspots);
        mHighConfidence = readUnnamedBed(mConfig.highConfidenceBed());

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.threads(), namedThreadFactory);
        mRefGenome = new IndexedFastaSequenceFile(new File(mConfig.refGenome()));
        mQualityRecalibrationSupplier = new QualityRecalibrationSupplier(mExecutorService, mRefGenome, mConfig);

        mVcfFile = new SageVCF(mRefGenome, mConfig);
        SG_LOGGER.info("Writing to file: {}", mConfig.outputFile());

        // Validate Coverage Bed
        if(mConfig.panelOnly() && !mCoveragePanel.isEmpty())
        {
            if(!GenomeRegionsValidation.isSubset(mPanelWithHotspots.values(), mCoveragePanel.values()))
            {
                throw new IOException("Coverage bed must be a subset of panel bed when running in panel only mode");
            }
        }
    }

    @NotNull
    private Coverage createCoverage()
    {
        Set<String> samples = Sets.newHashSet();
        if(!mConfig.coverageBed().isEmpty())
        {
            samples.addAll(mConfig.tumor());
        }
        return new Coverage(samples, mCoveragePanel.values());
    }

    private void run() throws InterruptedException, ExecutionException, IOException
    {
        long timeStamp = System.currentTimeMillis();
        final Coverage coverage = createCoverage();

        final Map<String, QualityRecalibrationMap> recalibrationMap = mQualityRecalibrationSupplier.get();
        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String contig = samSequenceRecord.getSequenceName();
            if(mConfig.chromosomes().isEmpty() || mConfig.chromosomes().contains(contig))
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
            String filename = mConfig.geneCoverageFile(sample);
            GeneDepthFile.write(filename, coverage.depth(sample));
        }

        long timeTaken = System.currentTimeMillis() - timeStamp;
        SG_LOGGER.info("Completed in {} seconds", timeTaken / 1000);
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = mConfig.referenceBam().isEmpty() ? mConfig.tumorBam().get(0) : mConfig.referenceBam().get(0);
        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.validationStringency())
                .referenceSource(new ReferenceSource(mRefGenome))
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
                mConfig,
                mExecutorService,
                mHotspots.get(chromosome),
                mPanelWithHotspots.get(chromosome),
                mHighConfidence.get(chromosome),
                qualityRecalibrationMap,
                coverage,
                mVcfFile::write);
    }

    @Override
    public void close() throws IOException
    {
        mVcfFile.close();
        mRefGenome.close();
        mExecutorService.shutdown();
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
        if(!mConfig.hotspots().isEmpty())
        {
            SG_LOGGER.info("Reading hotspot vcf: {}", mConfig.hotspots());
            return VariantHotspotFile.readFromVCF(mConfig.hotspots());
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
            SG_LOGGER.info("Reading bed file: {}", panelBed);
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
            SG_LOGGER.info("Reading bed file: {}", panelBed);
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