package com.hartwig.hmftools.svanalysis.gene;

import static com.hartwig.hmftools.svanalysis.types.SvaConfig.DATA_OUTPUT_DIR;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class GenerateEnsemblDataCache
{
    private static final Logger LOGGER = LogManager.getLogger(GenerateEnsemblDataCache.class);

    public static void writeEnsemblDataFiles(final CommandLine cmd)
    {
        EnsemblDAO ensemblData = new EnsemblDAO(cmd);
        ensemblData.writeDataCacheFiles(cmd.getOptionValue(DATA_OUTPUT_DIR));
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        Configurator.setRootLevel(Level.DEBUG);

        LOGGER.info("writing Ensembl data files");

        writeEnsemblDataFiles(cmd);

        LOGGER.info("Ensembl data cache complete");
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(DATA_OUTPUT_DIR, true, "Directory to write Ensembl data files");
        EnsemblDAO.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }


}