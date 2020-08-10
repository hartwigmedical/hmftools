package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.drivers.DriverAnnotation;
import com.hartwig.hmftools.cup.sample.SampleTraits;
import com.hartwig.hmftools.cup.sigs.SignatureAnnotation;
import com.hartwig.hmftools.cup.svs.SvAnnotation;

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
    private final SampleTraits mSampleTraits;
    private final SvAnnotation mSvAnnotation;

    private BufferedWriter mSampleDataWriter;

    public SampleAnalyser(final CommandLine cmd)
    {
        mConfig = new SampleAnalyserConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);

        mSnvSignatures = new SignatureAnnotation(mConfig, mSampleDataCache);
        mDrivers = new DriverAnnotation(mConfig, mSampleDataCache);
        mSampleTraits = new SampleTraits(mConfig, mSampleDataCache, mSnvSignatures);
        mSvAnnotation = new SvAnnotation(mConfig, mSampleDataCache);

        mSampleDataWriter = null;
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadSampleData(cmd.getOptionValue(SPECIFIC_SAMPLE_DATA), cmd.getOptionValue(SAMPLE_DATA_FILE));
        mSampleDataCache.loadReferenceSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));
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

            processSample(specificSample);
        }
        else
        {
            int sampleCount = 0;
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                processSample(sample);
                ++sampleCount;

                if((sampleCount % 100) == 0)
                {
                    CUP_LOGGER.info("processed {} samples", sampleCount);
                }
            }
        }

        closeBufferedWriter(mSampleDataWriter);

        CUP_LOGGER.info("CUP analysis complete");
    }

    private void initialiseOutputFiles()
    {
        try
        {
            final String outputFilename = mSampleDataCache.isSingleSample() ?
                    mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.data.csv"
                    : mConfig.formOutputFilename("SAMPLE_DATA");

            mSampleDataWriter = createBufferedWriter(outputFilename, false);

            mSampleDataWriter.write("SampleId,CancerType,Category,DataType,Value,RefCancerType,RefCancerTypeValue");
            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    private void processSample(final SampleData sampleData)
    {
        final List<SampleResult> traitsResults = mSampleTraits.processSample(sampleData);
        writeSampleData(sampleData, traitsResults);

        final List<SampleResult> snvResults = mSnvSignatures.processSample(sampleData);
        writeSampleData(sampleData, snvResults);

        final List<SampleResult> svResults = mSvAnnotation.processSample(sampleData);
        writeSampleData(sampleData, svResults);

        final List<SampleResult> driverResults = mDrivers.processSample(sampleData);
        writeSampleData(sampleData, driverResults);
    }

    private void writeSampleData(final SampleData sampleData, final List<SampleResult> results)
    {
        if(results.isEmpty() || mSampleDataWriter == null)
            return;

        try
        {
            for(SampleResult result : results)
            {
                final String sampleStr = String.format("%s,%s,%s,%s,%s",
                        sampleData.Id, sampleData.CancerType, result.Category, result.DataType, result.Value.toString());

                for(Map.Entry<String,Double> cancerValues : result.CancerTypeValues.entrySet())
                {
                    mSampleDataWriter.write(String.format("%s,%s,%.3g",
                            sampleStr, cancerValues.getKey(), cancerValues.getValue()));
                    mSampleDataWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample data: {}", e.toString());
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
