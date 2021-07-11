package com.hartwig.hmftools.sage.apps;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageQualityRecalibrationBedApplication implements AutoCloseable
{
    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String REF_GENOME = "ref_genome";
    private static final String OUT = "out";
    private static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;

    public static void main(final String[] args)
    {
        final Options options = createOptions();
        try(final SageQualityRecalibrationBedApplication application = new SageQualityRecalibrationBedApplication(options, args))
        {
            application.run();
        }
        catch(Exception e)
        {
            SG_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageQualityRecalibrationBedApplication", options);
            System.exit(1);
        }
    }

    private final IndexedFastaSequenceFile refGenome;
    private final int sampleSize;
    private final String out;

    private SageQualityRecalibrationBedApplication(final Options options, final String... args)
            throws ParseException, FileNotFoundException
    {
        final CommandLine cmd = createCommandLine(args, options);

        this.refGenome = new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME)));
        this.sampleSize = getConfigValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE);
        this.out = cmd.getOptionValue(OUT);
    }

    public void run() throws IOException
    {
        List<BaseRegion> regions = createRegions(refGenome, Sets.newHashSet(), sampleSize);
        List<GenomeRegion> genRegions = regions.stream().map(x -> GenomeRegions.create(x.Chromosome, x.start(), x.end())).collect(Collectors.toList());
        NamedBedFile.writeUnnamedBedFile(out, genRegions);
    }

    private static final int END_BUFFER = 1000000;
    private static final int REGION_SIZE = 100000;

    private static List<BaseRegion> createRegions(
            final IndexedFastaSequenceFile refGenome, final Set<String> chromosomes, int regionSubsetSize)
    {
        List<BaseRegion> result = Lists.newArrayList();

        for(final SAMSequenceRecord sequenceRecord : refGenome.getSequenceDictionary().getSequences())
        {
            final String chromosome = sequenceRecord.getSequenceName();

            if(!chromosomes.isEmpty() && !chromosomes.contains(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) || !HumanChromosome.fromString(chromosome).isAutosome())
                continue;

            int start = sequenceRecord.getSequenceLength() - END_BUFFER - regionSubsetSize;
            int end = sequenceRecord.getSequenceLength() - (END_BUFFER + 1);

            while(start < end)
            {
                result.add(new BaseRegion(chromosome, start, start + REGION_SIZE - 1));
                start += REGION_SIZE;
            }

        }
        return result;
    }


    @Override
    public void close() throws Exception
    {
        refGenome.close();
    }

    @NotNull
    static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(REF_GENOME, true, "Path to indexed ref genome fasta file");
        options.addOption(BQR_SAMPLE_SIZE, true, "BQR sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");

        options.addOption(OUT, true, "Path to output bed");

        return options;
    }

    @NotNull
    static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
