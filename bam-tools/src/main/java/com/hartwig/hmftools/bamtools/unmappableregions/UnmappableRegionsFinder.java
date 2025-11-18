package com.hartwig.hmftools.bamtools.unmappableregions;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMSequenceRecord;

public class UnmappableRegionsFinder
{
    private final UnmappableRegionsConfig mConfig;

    private final RefGenomeSource mRefGenomeSource;
    private final List<GenomeRegion> mUnmappableRegions = new ArrayList<>();

    private static final byte UNKNOWN_BASE = 'N';

    private static final String COL_CHROM = "chrom";
    private static final String COL_START = "start";
    private static final String COL_END = "end";

    public UnmappableRegionsFinder(UnmappableRegionsConfig config)
    {
        mConfig = config;

        mRefGenomeSource = RefGenomeSource.loadRefGenome(mConfig.RefGenome);
    }

    private GenomeRegion createUnmappableRegion(String chromosome, int start, int end)
    {
        BT_LOGGER.debug("Found unmappable region: {}:{}-{}", chromosome, start, end);
        return GenomeRegions.create(chromosome, start, end);
    }

    public void processChromosome(String chromosome)
    {
        BT_LOGGER.info("Processing chromosome: {}", chromosome);

        int currentUnmappedLength = 0;

        int chromosomeLength = mRefGenomeSource.getChromosomeLength(chromosome);
        for(int position = 1; position <= chromosomeLength; position++)
        {
            byte base = mRefGenomeSource.getBases(chromosome, position, position)[0];

            if(base == UNKNOWN_BASE)
            {
                currentUnmappedLength++;
            }
            else if(currentUnmappedLength > 0)
            {
                int unmappableStart = position - currentUnmappedLength;
                int unmappableEnd = position - 1;

                GenomeRegion unmappableRegion = createUnmappableRegion(chromosome, unmappableStart, unmappableEnd);
                mUnmappableRegions.add(unmappableRegion);

                currentUnmappedLength = 0;
            }
        }

        // Handle chromosome ending with Ns
        if(currentUnmappedLength > 0)
        {
            int unmappableEnd = chromosomeLength;
            int unmappableStart = unmappableEnd - currentUnmappedLength + 1;

            GenomeRegion unmappableRegion = createUnmappableRegion(chromosome, unmappableStart, unmappableEnd);
            mUnmappableRegions.add(unmappableRegion);
        }
    }

    public void writeRegions()
    {
        try
        {
            BT_LOGGER.info("Writing unmappable regions to file: {}", mConfig.OutputBedFile);

            BufferedWriter writer = createBufferedWriter(mConfig.OutputBedFile);
            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(COL_CHROM);
            header.add(COL_START);
            header.add(COL_END);

            writer.write(header.toString());
            writer.newLine();

            for(GenomeRegion region : mUnmappableRegions)
            {
                StringJoiner line = new StringJoiner(TSV_DELIM);

                line.add(region.chromosome());
                line.add(String.valueOf(region.start()));
                line.add(String.valueOf(region.end()));

                writer.write(line.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("Failed to write to file: {}", mConfig.OutputBedFile, e);
            System.exit(1);
        }
    }

    public void run()
    {
        BT_LOGGER.info("Finding unmappable regions");

        List<String> chromosomes = mRefGenomeSource.refGenomeFile()
                .getSequenceDictionary()
                .getSequences()
                .stream().map(SAMSequenceRecord::getSequenceName).toList();

        for(String chromosome : chromosomes)
        {
            boolean isHumanChromosome = HumanChromosome.contains(RefGenomeFunctions.stripChrPrefix(chromosome));

            if(!isHumanChromosome && !mConfig.IncludeAllChromosomes)
            {
                BT_LOGGER.trace("Skipping non-human chromosome: {}", chromosome);
                continue;
            }

            processChromosome(chromosome);
        }

        writeRegions();
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        UnmappableRegionsConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        UnmappableRegionsConfig config = new UnmappableRegionsConfig(configBuilder);
        new UnmappableRegionsFinder(config).run();
    }
}
