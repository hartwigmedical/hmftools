package com.hartwig.hmftools.pave.external;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.pave.Mappability.MAPPABILITY_BED;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class MappabilityBuilder
{
    private final String mOutputDir;
    private final String mInputFile;
    private final RefGenomeVersion mRefGenomeVersion;

    public MappabilityBuilder(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mInputFile = cmd.getOptionValue(MAPPABILITY_BED);
        mRefGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));
    }

    public void run()
    {
        if(mInputFile == null || !Files.exists(Paths.get(mInputFile)))
        {
            PV_LOGGER.error("missing input file, exiting");
            System.exit(1);
        }

        PV_LOGGER.info("reading mappability file({})", mInputFile);

        try
        {
            BufferedReader reader = createBufferedReader(mInputFile);

            String line = null;

            final Map<String,List<MapEntry>> chrEntries = Maps.newHashMap();

            List<MapEntry> entryList = null;
            String currentChromosome = "";
            int entryCount = 0;

            while((line = reader.readLine()) != null)
            {
                final String[] values = line.split("\t", -1); // eg: 1       0       10000   id-1    0.000000

                String chromosome = values[0];

                if(!HumanChromosome.contains(chromosome))
                    continue;

                if(!chromosome.equals(currentChromosome))
                {
                    if(entryList != null)
                    {
                        PV_LOGGER.debug("chromosome({}) has {} entries", currentChromosome, entryList.size());
                    }

                    currentChromosome = chromosome;
                    entryList = Lists.newArrayList();
                    chrEntries.put(chromosome, entryList);
                }

                entryList.add(new MapEntry(
                        new BaseRegion(Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2])),
                        Double.parseDouble(values[4])));

                ++entryCount;
            }

            PV_LOGGER.info("loaded {} mappability entries from file({})", entryCount, mInputFile);

            // write in chromosomal order
            String ouputFile = String.format("%s/mappability_150.%s.bed", mOutputDir, mRefGenomeVersion);
            BufferedWriter writer = createBufferedWriter(ouputFile, false);

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());
                List<MapEntry> entries = chrEntries.get(chrStr);

                for(MapEntry entry : entries)
                {
                    writer.write(String.format("%s\t%d\t%d\t%.4f", chrStr, entry.Region.start(), entry.Region.end(), entry.Mappability));
                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load mappability file({}): {}", mInputFile, e.toString());
        }

        PV_LOGGER.info("Mappability file build complete");
    }

    private class MapEntry
    {
        public final BaseRegion Region;
        public final double Mappability;

        public MapEntry(final BaseRegion region, final double mappability)
        {
            Region = region;
            Mappability = mappability;
        }

        public String toString() { return String.format("%s map(%.4f)", Region, Mappability); }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(MAPPABILITY_BED, true, "Mappability bed file");
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);
        setLogLevel(cmd);

        MappabilityBuilder mappabilityBuilder = new MappabilityBuilder(cmd);
        mappabilityBuilder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
