package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.StringJoiner;
import java.util.stream.Stream;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PeachApplication {
    static final String CHAIN_FILE_DELIM = " ";
    @NotNull
    private final PeachConfig config;
    public PeachApplication(@NotNull final PeachConfig config)
    {
        this.config = config;
    }

    public void run(){
        if(!config.isValid())
        {
            PCH_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }
        PCH_LOGGER.info("running PEACH");

        PCH_LOGGER.info("creating output directory");
        File outputDirectory = new File(config.outputDir);
        if (! outputDirectory.exists() && ! outputDirectory.mkdirs())
        {
            PCH_LOGGER.error("could not create output directory, exiting");
            System.exit(1);
        }

        String v38Vcf;
        if (config.doLiftOver)
        {
            PCH_LOGGER.info("create adjusted chain file");
            String adjustedChainFile = getAdjustedChainFile(config.chainFile);

            PCH_LOGGER.info("do lift over");
            v38Vcf = getExtendedFileName(config.vcfFile, "liftover", ".vcf");
            String rejectVcf = getExtendedFileName(config.vcfFile, "reject", ".vcf");
            ProcessBuilder pb = new ProcessBuilder(
                    "java",
                    "-jar",
                    config.picardJar,
                    "LiftoverVcf",
                    "CHAIN=" + adjustedChainFile,
                    "INPUT=" + config.vcfFile,
                    "OUTPUT=" + v38Vcf,
                    "REFERENCE_SEQUENCE=" + config.targetRefGenome,
                    "REJECT=" + rejectVcf,
                    "RECOVER_SWAPPED_REF_ALT=true",
                    "WRITE_ORIGINAL_POSITION=true",
                    "WRITE_ORIGINAL_ALLELES=true"
            );
            try
            {
                pb.inheritIO();
                Process process = pb.start();
                int exitCode = process.waitFor();
                PCH_LOGGER.info("Picard exit code: {}", exitCode);
            }
            catch(IOException e)
            {
                PCH_LOGGER.error("Picard LiftoverVcf failed: ");
                e.printStackTrace();
                System.exit(1);
            }
            catch(InterruptedException e)
            {
                PCH_LOGGER.error("Picard LiftoverVcf was interrupted");
                e.printStackTrace();
                System.exit(1);
            }
        }

        PCH_LOGGER.info("finished running PEACH");
    }

    private String getAdjustedChainFile(String chainFile) {
        String adjustedChainFile = getExtendedFileName(chainFile, "adjusted", ".over");
        try (
            Stream<String> lines = Files.lines(Paths.get(chainFile));
            PrintWriter pw = new PrintWriter(adjustedChainFile, StandardCharsets.UTF_8)
        )
        {
            lines.forEachOrdered(line-> pw.println(getAdjustedChainFileLine(line)));
        }
        catch (IOException e)
        {
            PCH_LOGGER.error("Could not create adjusted chain file: ");
            e.printStackTrace();
            System.exit(1);
        }

        return adjustedChainFile;
    }

    private String getAdjustedChainFileLine(String line) {
        if (line.startsWith("chain"))
        {
            String[] items = line.split(CHAIN_FILE_DELIM);
            StringJoiner newLineJoiner = new StringJoiner(CHAIN_FILE_DELIM);
            for (int i = 0; i < items.length; i++)
            {
                if (i == 2)
                    newLineJoiner.add(RefGenomeFunctions.stripChrPrefix(items[i]));
                else
                    newLineJoiner.add(items[i]);
            }
            return newLineJoiner.toString();
        }
        else
            return line;
    }

    private String getExtendedFileName(String originalFileName, String addition, String addBefore){
        String[] fileItems = originalFileName.split("/");
        String filename = fileItems[fileItems.length - 1];
        int extensionIndex = filename.indexOf(addBefore);
        return config.outputDir + filename.substring(0, extensionIndex) + "." + addition + filename.substring(extensionIndex);
    }

    public static void main(String[] args) {
        final Options options = PeachConfig.createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            PeachConfig config = new PeachConfig(cmd);
            PeachApplication peachApplication = new PeachApplication(config);
            peachApplication.run();
        }
        catch(ParseException e)
        {
            PCH_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PeachApplication", options);
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