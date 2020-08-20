package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupCalcs.addPercentileClassifier;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.feature.FeatureAnnotation;
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

    private final FeatureAnnotation mFeatures;
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
        mFeatures = new FeatureAnnotation(mConfig, mSampleDataCache);
        mSampleTraits = new SampleTraits(mConfig, mSampleDataCache);
        mSvAnnotation = new SvAnnotation(mConfig, mSampleDataCache);

        mSampleDataWriter = null;
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadSampleData(cmd.getOptionValue(SPECIFIC_SAMPLE_DATA), cmd.getOptionValue(SAMPLE_DATA_FILE));
        mSampleDataCache.loadReferenceSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));

        // mark any samples included in the ref data set so they can be excluded from self-comparison
        mSampleDataCache.SampleDataList.stream()
                .filter(x -> mSampleDataCache.RefSampleCancerTypeMap.containsKey(x.Id)).forEach(x -> x.setRefSample());
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        if(!checkAnnotators())
            return;

        initialiseOutputFiles();

        if(mSampleDataCache.SpecificSample != null)
        {
            final SampleData specificSample = mSampleDataCache.SpecificSample;

            CUP_LOGGER.info("sample({}) running CUP analysis", specificSample.Id);
            processSample(specificSample);
        }
        else
        {
            int sampleCount = 0;
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                // CUP_LOGGER.debug("sample({}) running CUP analysis", sample.Id);

                processSample(sample);

                if(!checkAnnotators())
                    break;

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

    private boolean checkAnnotators()
    {
        if(!mSvAnnotation.isValid() || !mFeatures.isValid() || !mSnvSignatures.isValid() || !mSampleTraits.isValid())
        {
            CUP_LOGGER.error("invalid init: traits({}) sigs({}) SVs({}) features({})",
                    mSampleTraits.isValid(), mSnvSignatures.isValid(), mSvAnnotation.isValid(), mFeatures.isValid());
            return false;
        }

        return true;
    }

    private void initialiseOutputFiles()
    {
        try
        {
            final String outputFilename = mSampleDataCache.isSingleSample() ?
                    mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.data.csv"
                    : mConfig.formOutputFilename("SAMPLE_DATA");

            mSampleDataWriter = createBufferedWriter(outputFilename, false);

            mSampleDataWriter.write("SampleId,CancerType,PrimaryLocation,CancerSubtype,Category,ResultType,DataType,Value,RefCancerType,RefValue");
            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    private void processSample(final SampleData sample)
    {
        final List<SampleResult> allResults = Lists.newArrayList();

        final List<SampleResult> traitsResults = mSampleTraits.processSample(sample);
        allResults.addAll(traitsResults);

        final List<SampleResult> snvResults = mSnvSignatures.processSample(sample);
        allResults.addAll(snvResults);

        final List<SampleResult> svResults = mSvAnnotation.processSample(sample);
        allResults.addAll(svResults);

        final List<SampleResult> driverResults = mFeatures.processSample(sample);
        allResults.addAll(driverResults);

        addPercentileClassifier(sample, allResults);

        writeSampleData(sample, allResults);
    }

    private void writeSampleData(final SampleData sampleData, final List<SampleResult> results)
    {
        if(results.isEmpty() || mSampleDataWriter == null)
            return;

        try
        {
            for(SampleResult result : results)
            {
                final String sampleStr = String.format("%s,%s,%s,%s,%s,%s,%s,%s",
                        sampleData.Id, sampleData.CancerType, sampleData.OriginalCancerType, sampleData.CancerSubtype,
                        result.Category, result.ResultType, result.DataType, result.Value.toString());

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
