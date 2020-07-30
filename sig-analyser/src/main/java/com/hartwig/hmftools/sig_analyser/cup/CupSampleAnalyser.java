package com.hartwig.hmftools.sig_analyser.cup;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CUP_LOGGER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_SIG_MIN_COUNT;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_SIG_MIN_PERCENT;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CupSampleAnalyser
{
    private final CupConfig mConfig;

    private SigMatrix mSampleCounts;
    private int[] mSampleTotals;

    private final List<CupSampleData> mSampleDataList;
    private final Map<String,List<CupSampleData>> mCancerSampleData;

    private final CupDrivers mDrivers;
    private final CupSnvSignatures mSnvSignatures;

    private BufferedWriter mSampleDataWriter;

    public CupSampleAnalyser(final CupConfig config)
    {
        mConfig = config;

        mSampleCounts = null;
        mSampleDataList = Lists.newArrayList();
        mCancerSampleData = Maps.newHashMap();
        mSampleTotals = null;

        mSnvSignatures = new CupSnvSignatures(config, mSampleDataList, mCancerSampleData);
        mDrivers = new CupDrivers(config, mSampleDataList);

        mSampleDataWriter = null;
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        final List<CupSampleData> sampleDataList = loadSampleData();

        if(!loadSampleCounts(sampleDataList))
            return;


        initialiseOutputFiles();

        mSnvSignatures.run();

        mDrivers.run();

        mSampleDataList.forEach(x -> writeSampleData(x));

        closeBufferedWriter(mSampleDataWriter);

        CUP_LOGGER.info("CUP SNV run complete");
    }


    private final List<CupSampleData> loadSampleData()
    {
        final List<CupSampleData> sampleDataList = Lists.newArrayList();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(mConfig.SampleDataFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileData)
            {
                final String[] items = line.split(",", -1);
                final String sampleId = items[fieldsIndexMap.get("SampleId")];
                final String cancerType = items[fieldsIndexMap.get("CancerType")];

                final String cancerSubtype = fieldsIndexMap.containsKey("CancerSubtype") ?
                        items[fieldsIndexMap.get("CancerSubtype")] : CANCER_SUBTYPE_OTHER;

                sampleDataList.add(new CupSampleData(sampleId, cancerType, cancerSubtype));
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample data file({}): {}", mConfig.SampleDataFile, e.toString());
        }

        return sampleDataList;
    }

    private boolean loadSampleCounts(final List<CupSampleData> sampleDataList)
    {
        final GenericDataCollection collection = GenericDataLoader.loadFile(mConfig.SnvSampleCountsFile);

        if(collection.getFieldNames().size() != sampleDataList.size())
        {
            CUP_LOGGER.error("sample count({} and data({}) mismatch");
            return false;
        }

        for(int s = 0; s < collection.getFieldNames().size(); ++s)
        {
            final String sampleId = collection.getFieldNames().get(s);
            CupSampleData sampleData = sampleDataList.stream().filter(x -> x.SampleId.equals(sampleId)).findFirst().orElse(null);

            if(sampleData == null)
                return false;

            sampleData.setSampleIndex(s);
            mSampleDataList.add(sampleData);

            List<CupSampleData> cancerSampleData = mCancerSampleData.get(sampleData.CancerType);
            if(cancerSampleData == null)
                mCancerSampleData.put(sampleData.CancerType, Lists.newArrayList(sampleData));
            else
                cancerSampleData.add(sampleData);
        }

        mSampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        mSampleCounts.cacheTranspose();

        return true;
    }

    private void initialiseOutputFiles()
    {
        try
        {
            mSampleDataWriter = createBufferedWriter(mConfig.formOutputFilename("css_data"), false);

            mSampleDataWriter.write("SampleId,CancerType,CancerSubtype");
            mSampleDataWriter.write(mSnvSignatures.getHeader());
            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    private void writeSampleData(final CupSampleData sampleData)
    {
        try
        {
            mSampleDataWriter.write(String.format("%s,%s,%s",
                    sampleData.SampleId, sampleData.CancerType, sampleData.CancerSubtype));

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
        CupConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        final CupConfig config = new CupConfig(cmd);

        CupSampleAnalyser sampleAnalyser = new CupSampleAnalyser(config);
        sampleAnalyser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
