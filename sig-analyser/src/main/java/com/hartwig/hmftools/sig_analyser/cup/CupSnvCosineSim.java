package com.hartwig.hmftools.sig_analyser.cup;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CUP_LOGGER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_CSS_THRESHOLD;

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
import com.hartwig.hmftools.sig_analyser.common.CosineSim;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.internal.$guava$.collect.$Interner;
import org.jetbrains.annotations.NotNull;

public class CupSnvCosineSim
{
    private final CupConfig mConfig;

    private SigMatrix mSampleCounts;
    private int[] mSampleTotals;

    private final List<CupSampleData> mSampleDataList;
    private final Map<String,List<CupSampleData>> mCancerSampleData;

    private BufferedWriter mCssPairsWriter;

    public CupSnvCosineSim(final CupConfig config)
    {
        mConfig = config;

        mSampleCounts = null;
        mSampleDataList = Lists.newArrayList();
        mCancerSampleData = Maps.newHashMap();
        mSampleTotals = null;
        mCssPairsWriter = null;
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

        int sampleCount = mSampleCounts.Cols;

        mSampleTotals = new int[sampleCount];

        for(int s = 0; s < sampleCount; ++s)
        {
            mSampleTotals[s] = (int)sumVector(mSampleCounts.getCol(s));
        }

        for(int s1 = 0; s1 < sampleCount; ++s1)
        {
            final CupSampleData sampleData1 = mSampleDataList.get(s1);
            final double[] sampleCounts1 = mSampleCounts.getCol(s1);

            // for(int s2 = s1 + 1; s2 < sampleCount; ++s2)
            for(int s2 = 0; s2 < sampleCount; ++s2)
            {
                if(s1 == s2)
                    continue;

                final CupSampleData sampleData2 = mSampleDataList.get(s2);
                final double[] sampleCounts2 = mSampleCounts.getCol(s2);

                double css = CosineSim.calcCSS(sampleCounts1, sampleCounts2);

                if(css < SNV_CSS_THRESHOLD)
                    continue;

                int cancerTypeCount = mCancerSampleData.get(sampleData2.CancerType).size();
                double weightedCss = pow((css - SNV_CSS_THRESHOLD) * 100, 2) / cancerTypeCount;

                sampleData1.addSampleCss(sampleData2.CancerType, weightedCss);

                // writeSampleCssPairData(mSampleIds.get(s1), mSampleIds.get(s2), css, mSampleTotals[s1], mSampleTotals[s2]);
            }

            writeSampleData(sampleData1);
        }

        closeBufferedWriter(mCssPairsWriter);

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
        final GenericDataCollection collection = GenericDataLoader.loadFile(mConfig.SnvSampleCounts);

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
            mCssPairsWriter = createBufferedWriter(mConfig.formOutputFilename("css_data"), false);

            mCssPairsWriter.write("SampleId,CancerType,CancerSubtype,SnvCount,TotalWeightedCss,MatchedCancerType,MatchedCancerCss");
            // mCssPairsWriter.write("SampleId1,SampleId2,CSS,SampleTotal1,SampleTotal2");
            mCssPairsWriter.newLine();
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
            final String sampleStr = String.format("%s,%s,%s,%d",
                    sampleData.SampleId, sampleData.CancerType, sampleData.CancerSubtype, mSampleTotals[sampleData.index()]);

            double totalWeightedCss = sampleData.getTotalWeightedCss();

            final Map<String,Double> cssDataMap = sampleData.getCancerCssTotals();

            if(cssDataMap.isEmpty())
            {
                mCssPairsWriter.write(String.format("%s,0,NONE,0", sampleStr));
                mCssPairsWriter.newLine();
                return;
            }

            for(Map.Entry<String,Double> cancerCssData : cssDataMap.entrySet())
            {
                mCssPairsWriter.write(String.format("%s,%.4f,%s,%.4f",
                        sampleStr, totalWeightedCss, cancerCssData.getKey(), cancerCssData.getValue()));
                mCssPairsWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }

    private void writeSampleCssPairData(final String sampleId1, final String sampleId2, double css, int sampleTotal1, int sampleTotal2)
    {
        try
        {
            mCssPairsWriter.write(String.format("%s,%s,%.6f,%d,%d", sampleId1, sampleId2, css, sampleTotal1, sampleTotal2));
            mCssPairsWriter.newLine();
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

        CupSnvCosineSim snvCosineSim = new CupSnvCosineSim(config);
        snvCosineSim.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
