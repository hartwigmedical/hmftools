package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleData;
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

    private final List<SampleData> mSampleDataList;
    private SampleData mSpecificSample;
    private final Map<String,List<SampleData>> mCancerSampleData;

    private final DriverAnnotation mDrivers;
    private final SignatureAnnotation mSnvSignatures;

    private BufferedWriter mSampleDataWriter;

    public SampleAnalyser(final CommandLine cmd)
    {
        mConfig = new SampleAnalyserConfig(cmd);

        mSampleDataList = Lists.newArrayList();
        mCancerSampleData = Maps.newHashMap();
        mSpecificSample = null;

        loadSampleData(cmd);

        mSnvSignatures = new SignatureAnnotation(mConfig, mSampleDataList, mCancerSampleData);
        mDrivers = new DriverAnnotation(mConfig, mSampleDataList);

        mSampleDataWriter = null;
    }

    private void loadSampleData(final CommandLine cmd)
    {
        if(cmd.hasOption(SPECIFIC_SAMPLE_DATA))
        {
            final String[] sampleItems = cmd.getOptionValue(SPECIFIC_SAMPLE_DATA).split(SUBSET_DELIM, -1);
            String sampleId = sampleItems[0];
            String cancerType = CANCER_TYPE_UNKNOWN;
            String cancerSubtype = CANCER_SUBTYPE_OTHER;

            if(sampleItems.length >= 2)
                cancerType = sampleItems[1];

            if(sampleItems.length == 3)
                cancerSubtype = sampleItems[2];

            mSpecificSample = new SampleData(sampleId, cancerType, cancerSubtype);
            mSampleDataList.add(mSpecificSample);
        }

        if(cmd.hasOption(SAMPLE_DATA_FILE))
        {
            try
            {
                final List<String> fileData = Files.readAllLines(new File(mConfig.SampleDataFile).toPath());

                final String header = fileData.get(0);
                fileData.remove(0);

                final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

                for(final String line : fileData)
                {
                    final String[] items = line.split(DATA_DELIM, -1);
                    final String sampleId = items[fieldsIndexMap.get("SampleId")];
                    final String cancerType = items[fieldsIndexMap.get("CancerType")];

                    final String cancerSubtype = fieldsIndexMap.containsKey("CancerSubtype") ?
                            items[fieldsIndexMap.get("CancerSubtype")] : CANCER_SUBTYPE_OTHER;

                    if(mSpecificSample != null && mSpecificSample.Id.equals(sampleId))
                        continue;

                    mSampleDataList.add(new SampleData(sampleId, cancerType, cancerSubtype));
                }
            }
            catch (IOException e)
            {
                CUP_LOGGER.error("failed to read sample data file({}): {}", mConfig.SampleDataFile, e.toString());
            }
        }

        // build a cache of samples per cancer type
        for(SampleData sampleData : mSampleDataList)
        {
            List<SampleData> cancerSampleData = mCancerSampleData.get(sampleData.CancerType);
            if(cancerSampleData == null)
                mCancerSampleData.put(sampleData.CancerType, Lists.newArrayList(sampleData));
            else
                cancerSampleData.add(sampleData);
        }
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        initialiseOutputFiles();

        if(mSpecificSample != null)
        {
            mSnvSignatures.processSample(mSpecificSample);
            mDrivers.processSample(mSpecificSample);
            writeSampleData(mSpecificSample);
        }
        else
        {
            for(SampleData sample : mSampleDataList)
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
