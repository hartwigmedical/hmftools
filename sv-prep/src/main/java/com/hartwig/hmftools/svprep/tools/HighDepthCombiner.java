package com.hartwig.hmftools.svprep.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HighDepthCombiner
{
    private final List<String> mInputFiles;
    private final BufferedWriter mWriter;

    private static final String INPUT_FILES = "input_files";
    private static final String OUTPUT_FILE = "output_file";

    public HighDepthCombiner(final CommandLine cmd)
    {
        mInputFiles = Arrays.stream(cmd.getOptionValue(INPUT_FILES).split(",")).collect(Collectors.toList());
        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_FILE));
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            SV_LOGGER.error("no input files specified");
            System.exit(1);
        }


        closeBufferedWriter(mWriter);

        SV_LOGGER.info("high depth discovery complete");

        // write output VCF
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        SV_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,PosStart,PosEnd,BaseDepth");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to initialise writer: {}", e.toString());
        }

        return null;
    }

    public void writeCombinedResults(final BufferedWriter writer, final List<HighDepthRegion> regions)
    {
        try
        {
            for(HighDepthRegion region : regions)
            {
                writer.write(format("%s,%d,%d,%d", region.Region.Chromosome, region.Region.start(), region.Region.end(), region.Depth));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write region: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(cmd);
        highDepthCombiner.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
