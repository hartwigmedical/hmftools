package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.LOG_DEBUG;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.RUN_MODE_BOTH;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.RUN_MODE_POST_PROCESS;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.RUN_MODE_VCF_PARSE;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.createCommandLine;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.createOptions;

import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class BachelorApplication {

    private final GermlineVcfParser mGermlineVcfParser;
    private final BachelorPostProcess mPostProcessor;

    private final BachelorConfig mConfig;

    private final List<RunDirectory> mSampleDataDirectories;


    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);

    private BachelorApplication(final CommandLine cmd)
    {
        mSampleDataDirectories = Lists.newArrayList();

        mConfig = new BachelorConfig(cmd);

        setSampleDataDirectories();

        if(mConfig.RunMode.equals(RUN_MODE_BOTH) || mConfig.RunMode.equals(RUN_MODE_VCF_PARSE))
        {
            mGermlineVcfParser = new GermlineVcfParser(mConfig, cmd);
        }
        else
        {
            mGermlineVcfParser = null;
        }

        if(mConfig.RunMode.equals(RUN_MODE_BOTH) || mConfig.RunMode.equals(RUN_MODE_POST_PROCESS))
        {
            mPostProcessor = new BachelorPostProcess(mConfig, cmd);
        }
        else
        {
            mPostProcessor = null;
        }
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            LOGGER.error("exiting due to invalid config");
            return;
        }

        for (int i = 0; i < mSampleDataDirectories.size(); ++i)
        {
            final RunDirectory runDir = mSampleDataDirectories.get(i);

            String sampleId = "";

            if(!mConfig.IsBatchMode)
            {
                sampleId = mConfig.SampleId;
            }
            else
            {
                try
                {
                    final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDir.sampleDir().toString());
                    sampleId = runContext.tumorSample();
                }
                catch (Exception e)
                {
                    // Skip using meta data
                    // sampleId = runDir.getPatientID();
                    continue;
                }
            }

            if(!mConfig.RestrictedSampleIds.isEmpty() && !mConfig.RestrictedSampleIds.contains(sampleId))
            {
                LOGGER.info("Skipping sampleId({}) not in specified list", sampleId);
                continue;
            }

            if(mGermlineVcfParser != null)
            {
                mGermlineVcfParser.run(runDir, sampleId, mConfig.OutputDir);
            }

            if(mPostProcessor != null)
            {
                List<BachelorGermlineVariant> bachelorRecords = mGermlineVcfParser != null ? mGermlineVcfParser.getBachelorRecords() : null;
                mPostProcessor.run(bachelorRecords);
            }

            if(mConfig.MaxBatchDirectories > 0 && i >= mConfig.MaxBatchDirectories)
                break;
        }

        if(mGermlineVcfParser != null)
        {
            mGermlineVcfParser.close();
        }

        if(mPostProcessor != null)
        {
            mPostProcessor.close();
        }

        LOGGER.info("bachelor run complete");
    }

    private void setSampleDataDirectories()
    {
        final Path sampleDataPath = Paths.get(mConfig.SampleDataDir);

        if(mConfig.IsBatchMode)
        {
            try (final Stream<Path> stream = Files.walk(sampleDataPath, 1, FileVisitOption.FOLLOW_LINKS).parallel())
            {
                mSampleDataDirectories.addAll(stream.filter(p -> p.toFile().isDirectory())
                        .filter(p -> !p.equals(sampleDataPath))
                        .map(RunDirectory::new)
                        .collect(Collectors.toList()));
            }
            catch (Exception e)
            {
                LOGGER.error("failed find directories for batch run: {}", e.toString());
            }

            LOGGER.info("Found {} batch directories", mSampleDataDirectories.size());
        }
        else
        {
            mSampleDataDirectories.add(new RunDirectory(sampleDataPath));
        }
    }

    public static void main(final String... args)
    {
        final Options options = createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            BachelorApplication bachelorApp = new BachelorApplication(cmd);
            bachelorApp.run();
        }
        catch (final ParseException e)
        {
            printHelpAndExit(options);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    private static void printHelpAndExit(final Options options)
    {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bachelor", "Determines eligibility", options, "", true);
        System.exit(1);
    }
}
