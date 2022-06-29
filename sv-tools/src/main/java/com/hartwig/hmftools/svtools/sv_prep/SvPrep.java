package com.hartwig.hmftools.svtools.sv_prep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svtools.sv_prep.SvConfig.createCmdLineOptions;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SvPrep
{
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;


    public SvPrep(final CommandLine cmd)
    {
        mConfig = new SvConfig(cmd);
        mWriter = new ResultsWriter(mConfig);
    }

    public void run()
    {
        SV_LOGGER.info("sample({}) starting SV prep", mConfig.SampleId);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            SV_LOGGER.info("processing chromosome({})", chromosomeStr);

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosomeStr, mConfig, mWriter);
            chromosomeTask.process();
            System.gc();
        }

        mWriter.close();

        SV_LOGGER.info("SV prep complete");
    }

    public static void main(@NotNull final String[] args) throws Exception
    {
        // final VersionInfo version = new VersionInfo("??.version");
        // SV_LOGGER.info("?? version: {}", version.version());

        final Options options = createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            SvPrep svPrep = new SvPrep(cmd);
            svPrep.run();
        }
        catch(ParseException e)
        {
            SV_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Isofox", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
