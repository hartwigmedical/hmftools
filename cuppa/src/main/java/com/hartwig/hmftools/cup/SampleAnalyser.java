package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.drivers.DriverAnnotation;
import com.hartwig.hmftools.cup.sigs.SignatureAnnotation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SampleAnalyser
{
    private final SampleAnalyserConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final DriverAnnotation mDrivers;
    private final SignatureAnnotation mSnvSignatures;

    private BufferedWriter mSampleDataWriter;

    public SampleAnalyser(final CommandLine cmd)
    {
        mConfig = new SampleAnalyserConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);

        mSnvSignatures = new SignatureAnnotation(mConfig, mSampleDataCache);
        mDrivers = new DriverAnnotation(mConfig, mSampleDataCache);

        mSampleDataWriter = null;
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadSampleData(cmd.getOptionValue(SPECIFIC_SAMPLE_DATA), cmd.getOptionValue(SAMPLE_DATA_FILE));
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        initialiseOutputFiles();

        if(mSampleDataCache.SpecificSample != null)
        {
            final SampleData specificSample = mSampleDataCache.SpecificSample;
            mSnvSignatures.processSample(specificSample);
            mDrivers.processSample(specificSample);
            writeSampleData(specificSample);
        }
        else
        {
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                mSnvSignatures.processSample(sample);

                mDrivers.processSample(sample);

                writeSampleData(sample);
            }
        }

        mSnvSignatures.close();
        mDrivers.close();

        closeBufferedWriter(mSampleDataWriter);

        CUP_LOGGER.info("CUP analysis complete");
    }

    private void initialiseOutputFiles()
    {
        try
        {
            mSampleDataWriter = createBufferedWriter(mConfig.formOutputFilename("SAMPLE_DATA"), false);

            mSampleDataWriter.write("SampleId,CancerType,CancerSubtype");
            mSampleDataWriter.write("," + mSnvSignatures.getHeader());
            mSampleDataWriter.write("," + mDrivers.getHeader());
            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    private void writeSampleData(final SampleData sampleData)
    {
        try
        {
            mSampleDataWriter.write(String.format("%s,%s,%s",
                    sampleData.Id, sampleData.CancerType, sampleData.CancerSubtype));

            mSampleDataWriter.write(String.format(",%s,%s",
                    mSnvSignatures.getSampleOutput(sampleData), mDrivers.getSampleOutput(sampleData)));

            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        SampleAnalyserConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        SampleAnalyser sampleAnalyser = new SampleAnalyser(cmd);
        sampleAnalyser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
