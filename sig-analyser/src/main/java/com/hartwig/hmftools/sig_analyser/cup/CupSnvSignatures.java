package com.hartwig.hmftools.sig_analyser.cup;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CUP_LOGGER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_SIG_MIN_COUNT;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.SNV_SIG_MIN_PERCENT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.sig_analyser.common.CosineSim;

public class CupSnvSignatures
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
    private BufferedWriter mSampleDataWriter;

    public CupSnvSignatures(final CupConfig config, final List<CupSampleData> sampleDataList, final Map<String,List<CupSampleData>> cancerSampleData)
    {
        mConfig = config;

        mSampleCounts = null;
        mSampleDataList = sampleDataList;
        mCancerSampleData = cancerSampleData;
        mSampleTotals = null;

        mSnvSignatures = null;
        mLeastSquaresFit = null;
        mSnvSigNames = Lists.newArrayList();

        mSampleDataWriter = null;
        mCssPairsWriter = null;

        loadSnvSignatures();
        initialiseOutputFile();
    }

    public void run()
    {
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
        }

        mSampleDataList.forEach(x -> writeSampleData(x));
        closeBufferedWriter(mSampleDataWriter);
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

    public String getHeader()
    {
        return "SnvCount,TotalWeightedCss,TopMatchCancerType,TopMatchCancerCss,SigData";
    }

    public String getSampleOutput(final CupSampleData sampleData)
    {
        String sampleDataStr = String.format("%d", mSampleTotals[sampleData.index()]);

        double totalWeightedCss = sampleData.getTotalWeightedCss();

        String sigDataStr = "";
        for(Map.Entry<String,Double> entry : sampleData.getSnvSigAllocations().entrySet())
        {
            sigDataStr = appendStr(sigDataStr, String.format("%s=%.2f", entry.getKey(), entry.getValue()), ';');
        }

        final Map<String,Double> cssDataMap = sampleData.getCancerCssTotals();

        String cssDataStr;
        if(cssDataMap.isEmpty())
        {
            cssDataStr = "0,NONE,0";
        }
        else
        {
            String topType = "";
            double topCss = 0;

            for(Map.Entry<String,Double> cancerCssData : cssDataMap.entrySet())
            {
                if(cancerCssData.getValue() > topCss)
                {
                    topCss = cancerCssData.getValue();
                    topType = cancerCssData.getKey();
                }
            }

            cssDataStr = String.format("%.4f,%s,%.4f", totalWeightedCss, topType, topCss);
        }

        return String.format("%s,%s,%s", sampleDataStr,cssDataStr,sigDataStr);
    }

    private void initialiseOutputFile()
    {
        try
        {
            mSampleDataWriter = createBufferedWriter(mConfig.formOutputFilename("SNV"), false);

            mSampleDataWriter.write("SampleId,SnvCount,TotalWeightedCss,MatchedCancerType,MatchedCancerCss,SigData");
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
            final String sampleStr = String.format("%s,%d",
                    sampleData.SampleId, mSampleTotals[sampleData.index()]);

            double totalWeightedCss = sampleData.getTotalWeightedCss();

            String sigDataStr = "";
            for(Map.Entry<String,Double> entry : sampleData.getSnvSigAllocations().entrySet())
            {
                sigDataStr = appendStr(sigDataStr, String.format("%s=%.2f", entry.getKey(), entry.getValue()), ';');
            }

            final Map<String,Double> cssDataMap = sampleData.getCancerCssTotals();

            if(cssDataMap.isEmpty())
            {
                mSampleDataWriter.write(String.format("%s,0,NONE,0,%s", sampleStr, sigDataStr));
                mSampleDataWriter.newLine();
                return;
            }

            for(Map.Entry<String,Double> cancerCssData : cssDataMap.entrySet())
            {
                mSampleDataWriter.write(String.format("%s,%.4f,%s,%.4f,%s",
                        sampleStr, totalWeightedCss, cancerCssData.getKey(), cancerCssData.getValue(), sigDataStr));
                mSampleDataWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
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
}
