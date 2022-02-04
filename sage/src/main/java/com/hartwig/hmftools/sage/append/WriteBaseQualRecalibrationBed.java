package com.hartwig.hmftools.sage.append;

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
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class WriteBaseQualRecalibrationBed
{
    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String REF_GENOME = "ref_genome";
    private static final String OUT = "out";
    private static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;

    private final IndexedFastaSequenceFile mRefGenome;
    private final int mSampleSize;
    private final String mOutputFile;

    private WriteBaseQualRecalibrationBed(final Options options, final String... args)
            throws ParseException, FileNotFoundException
    {
        final CommandLine cmd = createCommandLine(args, options);

        mRefGenome = new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME)));
        mSampleSize = getConfigValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE);
        mOutputFile = cmd.getOptionValue(OUT);
    }

    public void run() throws IOException
    {
        List<ChrBaseRegion> regions = createRegions(mRefGenome, Sets.newHashSet(), mSampleSize);

        SG_LOGGER.info("writing {} BQR regions to BED file: {}", regions.size(), mOutputFile);

        List<GenomeRegion> genRegions = regions.stream().map(x -> GenomeRegions.create(x.Chromosome, x.start(), x.end())).collect(Collectors.toList());
        NamedBedFile.writeUnnamedBedFile(mOutputFile, genRegions);

        SG_LOGGER.info("BED file write complete");

        mRefGenome.close();
    }

    private static final int END_BUFFER = 1000000;
    private static final int REGION_SIZE = 100000;

    private static List<ChrBaseRegion> createRegions(
            final IndexedFastaSequenceFile refGenome, final Set<String> chromosomes, int regionSubsetSize)
    {
        List<ChrBaseRegion> result = Lists.newArrayList();

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
                result.add(new ChrBaseRegion(chromosome, start, start + REGION_SIZE - 1));
                start += REGION_SIZE;
            }

        }
        return result;
    }

    public static void main(final String[] args)
    {
        final Options options = createOptions();
        try
        {
            final WriteBaseQualRecalibrationBed application = new WriteBaseQualRecalibrationBed(options, args);
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
