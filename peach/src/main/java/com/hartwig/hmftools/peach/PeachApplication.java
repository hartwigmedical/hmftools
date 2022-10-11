package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PeachApplication {
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

        String v38Vcf;
        if (config.doLiftOver)
        {
            PCH_LOGGER.info("do lift over");
            v38Vcf = getExtendedVcfFileName(config.vcfFile, "38");
            String rejectVcf = getExtendedVcfFileName(config.vcfFile, "reject");

            ProcessBuilder pb = new ProcessBuilder(
                    "java",
                    "-jar",
                    config.picardJar,
                    "--CHAIN",
                    config.chainFile,
                    "--INPUT",
                    config.vcfFile,
                    "--OUTPUT",
                    v38Vcf,
                    "--REFERENCE_SEQUENCE",
                    config.targetRefGenome,
                    "--REJECT",
                    rejectVcf
            );
            try {
                Process process = pb.start();
                int exitCode = process.waitFor();
                PCH_LOGGER.info("Picard exit code: {}", exitCode);
            }
            catch(IOException e){
                PCH_LOGGER.error("Picard LiftoverVcf failed: ");
                e.printStackTrace();
                System.exit(1);
            }
            catch(InterruptedException e){
                PCH_LOGGER.error("Picard LiftoverVcf was interrupted");
                System.exit(1);
            }
        }

        PCH_LOGGER.info("finished running PEACH");
    }

    private String getExtendedVcfFileName(String originalFileName, String extension){
        String[] fileItems = originalFileName.split("/");
        String filename = fileItems[fileItems.length - 1];
        int extensionIndex = filename.indexOf(".vcf");
        return config.outputDir + filename.substring(0, extensionIndex) + "." + extension + filename.substring(extensionIndex);
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