package com.hartwig.hmftools.virusdetect;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.virusdetect.VirusConstants.APP_NAME;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VirusApplication
{
    private final VirusConfig mConfig;
    private final ViralReference mViralReference;

    private static final Logger LOGGER = LogManager.getLogger(VirusApplication.class);

    public VirusApplication(VirusConfig config)
    {
        mConfig = config;

        LOGGER.info("Loading prerequisite data");
        mViralReference = ViralReference.load(config.viralRefFile(), config.viralRefInfoFile());
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting VirusDetect");
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        // TODO: placeholder pipeline; each step is replaced by its implementation as it lands.
        LOGGER.info("Extracting candidate viral reads -> FASTA (stub)");
        LOGGER.info("Realigning candidates to viral reference -> BAM (stub)");
        LOGGER.info("Computing per-contig stats over aligned BAM -> debug TSV (stub)");
        LOGGER.info("Selecting representative contig per oncology group (stub)");
        LOGGER.info("Filtering aligned BAM to representatives -> BAM (stub)");
        LOGGER.info("Computing per-contig stats over representative BAM (stub)");
        LOGGER.info("Annotating QC and writing detected TSV (stub)");

        LOGGER.info("VirusDetect complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        VirusConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try
        {
            VirusConfig config = VirusConfig.fromConfigBuilder(configBuilder);
            VirusApplication virusDetect = new VirusApplication(config);
            virusDetect.run();
        }
        catch(UserInputError e)
        {
            LOGGER.error("Bad input data: {}", e.getMessage());
            exit(1);
        }
        catch(IOException e)
        {
            LOGGER.error("IO error", e);
            exit(1);
        }
        catch(RuntimeException e)
        {
            LOGGER.error("Runtime error", e);
            exit(1);
        }
    }
}
