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
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
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
import org.jetbrains.annotations.NotNull;

public class CupSnvCosineSim
{
    private final CupConfig mConfig;

    private SigMatrix mSampleCounts;
    private int[] mSampleTotals;

    private final List<CupSampleData> mSampleDataList;
    private final Map<String,List<CupSampleData>> mCancerSampleData;

    private SigMatrix mSnvSignatures;
    private final List<String> mSnvSigNames;
    private LeastSquaresFit mLeastSquaresFit;

    private BufferedWriter mCssPairsWriter;

    public CupSnvCosineSim(final CupConfig config)
    {
        mConfig = config;

        mSampleCounts = null;
        mSampleDataList = Lists.newArrayList();
        mCancerSampleData = Maps.newHashMap();
        mSampleTotals = null;

        mSnvSignatures = null;
        mLeastSquaresFit = null;
        mSnvSigNames = Lists.newArrayList();

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

        loadSnvSignatures();

        initialiseOutputFiles();

        int sampleCount = mSampleCounts.Cols;

        mSampleTotals = new int[sampleCount];

        for(int s = 0; s < sampleCount; ++s)
        {
            mSampleTotals[s] = (int)sumVector(mSampleCounts.getCol(s));
        }

        for(int s = 0; s < sampleCount; ++s)
        {
            final CupSampleData sampleData = mSampleDataList.get(s);
            final double[] sampleCounts = mSampleCounts.getCol(s);

            // for(int s2 = s1 + 1; s2 < sampleCount; ++s2)
            for(int s2 = 0; s2 < sampleCount; ++s2)
            {
                if(s == s2)
                    continue;

                final CupSampleData otherSampleData = mSampleDataList.get(s2);

                if(otherSampleData.isUnknownCancerType())
                    continue;

                final double[] otherSampleCounts = mSampleCounts.getCol(s2);

                double css = CosineSim.calcCSS(sampleCounts, otherSampleCounts);

                if(css < SNV_CSS_THRESHOLD)
                    continue;

                int cancerTypeCount = mCancerSampleData.get(otherSampleData.CancerType).size();
                double weightedCss = pow((css - SNV_CSS_THRESHOLD) * 100, 2) / cancerTypeCount;

                sampleData.addSampleCss(otherSampleData.CancerType, weightedCss);

                // writeSampleCssPairData(mSampleIds.get(s1), mSampleIds.get(s2), css, mSampleTotals[s1], mSampleTotals[s2]);
            }

            fitSnvSignatures(sampleData);

            writeSampleData(sampleData);
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

    private void fitSnvSignatures(final CupSampleData sampleData)
    {
        if(mSnvSignatures == null)
            return;

        final double[] sampleCounts = mSampleCounts.getCol(sampleData.index());

        mLeastSquaresFit.initialise(mSnvSignatures.getData(), sampleCounts);
        mLeastSquaresFit.solve();

        final double[] sigAllocs = mLeastSquaresFit.getContribs();

        double sampleTotal = mSampleTotals[sampleData.index()];

        for(int sig = 0; sig < sigAllocs.length; ++sig)
        {
            double sigAlloc = sigAllocs[sig];
            double allocPerc = sigAlloc / sampleTotal;

            if(sigAlloc >= SNV_SIG_MIN_COUNT && allocPerc >= SNV_SIG_MIN_PERCENT)
            {
                sampleData.getSnvSigAllocations().put(mSnvSigNames.get(sig), allocPerc);
            }
        }
    }

    private void loadSnvSignatures()
    {
        if(mConfig.SnvSignaturesFile.isEmpty())
            return;

        // cosmic_sig_subset.csv
        final GenericDataCollection sigsCollection = GenericDataLoader.loadFile(mConfig.SnvSignaturesFile);
        mSnvSigNames.addAll(sigsCollection.getFieldNames());
        mSnvSignatures = DataUtils.createMatrixFromListData(sigsCollection.getData());
        mLeastSquaresFit = new LeastSquaresFit(mSnvSignatures.Rows, mSnvSignatures.Cols);
    }

    private void initialiseOutputFiles()
    {
        try
        {
            mCssPairsWriter = createBufferedWriter(mConfig.formOutputFilename("css_data"), false);

            mCssPairsWriter.write("SampleId,CancerType,CancerSubtype,SnvCount,TotalWeightedCss,MatchedCancerType,MatchedCancerCss,SigData");
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

            String sigDataStr = "";
            for(Map.Entry<String,Double> entry : sampleData.getSnvSigAllocations().entrySet())
            {
                sigDataStr = appendStr(sigDataStr, String.format("%s=%.2f", entry.getKey(), entry.getValue()), ';');
            }

            final Map<String,Double> cssDataMap = sampleData.getCancerCssTotals();

            if(cssDataMap.isEmpty())
            {
                mCssPairsWriter.write(String.format("%s,0,NONE,0,%s", sampleStr, sigDataStr));
                mCssPairsWriter.newLine();
                return;
            }

            for(Map.Entry<String,Double> cancerCssData : cssDataMap.entrySet())
            {
                mCssPairsWriter.write(String.format("%s,%.4f,%s,%.4f,%s",
                        sampleStr, totalWeightedCss, cancerCssData.getKey(), cancerCssData.getValue(), sigDataStr));
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
