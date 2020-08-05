package com.hartwig.hmftools.cup.sigs;

import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV_SIG;
import static com.hartwig.hmftools.cup.common.CupConstants.SNV_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.sigs.RefSignatures.populateRefSigContributions;

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
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.sig_analyser.common.CosineSim;

public class SignatureAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private SigMatrix mRefSampleCounts;
    private final List<String> mRefSampleNames;
    private final Map<String,Map<String,double[]>> mRefCancerSigContribPercentiles;

    private SigMatrix mSampleCounts;
    private final Map<String,Integer> mSampleCountsIndex;
    private final Map<String,Map<String,Double>> mSampleSigContributions;

    private BufferedWriter mCssPairsWriter;
    private BufferedWriter mSampleDataWriter;

    public SignatureAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleCounts = null;
        mSampleSigContributions = Maps.newHashMap();
        mSampleCountsIndex = Maps.newHashMap();

        mRefSampleCounts = null;
        mRefSampleNames = Lists.newArrayList();
        mRefCancerSigContribPercentiles = Maps.newHashMap();

        mSampleDataWriter = null;
        mCssPairsWriter = null;

        loadRefSigContribPercentiles();
        loadRefSampleCounts();
        loadSampleCounts();
        loadSigContributions();
        initialiseOutputFile();
    }

    private boolean loadRefSampleCounts()
    {
        final GenericDataCollection collection = GenericDataLoader.loadFile(mConfig.RefSnvCountsFile);

        mRefSampleNames.addAll(collection.getFieldNames());
        mRefSampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        mRefSampleCounts.cacheTranspose();

        return true;
    }
    private boolean loadSampleCounts()
    {
        final GenericDataCollection collection = GenericDataLoader.loadFile(mConfig.SampleSnvCountsFile);

        for(int s = 0; s < collection.getFieldNames().size(); ++s)
        {
            final String sampleId = collection.getFieldNames().get(s);
            mSampleCountsIndex.put(sampleId, s);
        }

        mSampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        mSampleCounts.cacheTranspose();

        return true;
    }

    private void loadSigContributions()
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(mConfig.SampleSigContribFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String sampleId = items[0];
                String sigName = items[1];
                double sigContrib = Double.parseDouble(items[2]);

                Map<String,Double> sigContribs = mSampleSigContributions.get(sampleId);

                if(sigContribs == null)
                {
                    sigContribs = Maps.newHashMap();
                    mSampleSigContributions.put(sampleId, sigContribs);
                }

                sigContribs.put(sigName, sigContrib);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contribution data file({}): {}", mConfig.SampleSigContribFile, e.toString());
        }
    }

    public void processCohort()
    {
        for(SampleData sample : mSampleDataCache.SampleDataList)
        {
            processSample(sample);
        }

        closeBufferedWriter(mSampleDataWriter);
    }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = org.apache.commons.compress.utils.Lists.newArrayList();

        Integer sampleCountsIndex = mSampleCountsIndex.get(sample.Id);

        if(sampleCountsIndex == null)
        {
            CUP_LOGGER.error("sample({}) has no SNV data", sample.Id);
            return results;
        }

        addCssResults(sample, sampleCountsIndex, results);

        addSigContributionResults(sample, sampleCountsIndex, results);

        // fitSnvSignatures(sample, sampleCountsIndex);

        writeSampleData(sample);

        return results;
    }

    private void addCssResults(final SampleData sample, int sampleCountsIndex, final List<SampleResult> results)
    {
        final double[] sampleCounts = mSampleCounts.getCol(sampleCountsIndex);

        sample.setSnvCount((int)sumVector(sampleCounts));

        int refSampleCount = mRefSampleCounts.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        for(int s = 0; s < refSampleCount; ++s)
        {
            final String refSampleId = mRefSampleNames.get(s);

            if(refSampleId.equals(sample.Id))
                continue;

            final double[] otherSampleCounts = mRefSampleCounts.getCol(s);

            double css = CosineSim.calcCSS(sampleCounts, otherSampleCounts);

            if(css < SNV_CSS_THRESHOLD)
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            int cancerTypeCount = mSampleDataCache.RefCancerSampleData.get(refCancerType).size();
            double weightedCss = pow((css - SNV_CSS_THRESHOLD) * 100, 2) / cancerTypeCount;

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);

            // writeSampleCssPairData(mSampleIds.get(s1), mSampleIds.get(s2), css, mSampleTotals[s1], mSampleTotals[s2]);
        }

        results.add(new SampleResult(sample.Id, SNV, "CSS", 0, cancerCssTotals));
    }

    private void addSigContributionResults(final SampleData sample, int sampleCountsIndex, final List<SampleResult> results)
    {
        final Map<String,Double> sampleSigContribs = mSampleSigContributions.get(sample.Id);

        if(sampleSigContribs == null)
        {
            CUP_LOGGER.error("sample({}) sig contributions not found", sample.Id);
            return;
        }

        for(Map.Entry<String,Double> entry : sampleSigContribs.entrySet())
        {
            final String sigName = entry.getKey();
            double sampleSigContrib = entry.getValue();

            if(sampleSigContrib < 1)
                continue;

            for(Map.Entry<String,Map<String,double[]>> cancerContribs : mRefCancerSigContribPercentiles.entrySet())
            {
                final String cancerType = cancerContribs.getKey();
                final double[] refSigPercentiles = cancerContribs.getValue().get(sigName);

                if(refSigPercentiles == null)
                {
                    CUP_LOGGER.error("missing sig({}) data for cancerType({})", sigName, cancerType);
                    return;
                }

                Map<String,Double> cancerResults = Maps.newHashMap();
                double percentile = getPercentile(refSigPercentiles, sampleSigContrib);
                cancerResults.put(cancerType, percentile);

                results.add(new SampleResult(sample.Id, SNV_SIG, sigName.toUpperCase(), round(sampleSigContrib), cancerResults));
            }
        }
    }

    private String formatSigContribResult(double contrib)
    {
        return String.format("%.1f", contrib);
    }

    public String getHeader()
    {
        return "SnvCount,TotalWeightedCss,TopMatchCancerType,TopMatchCancerCss,SigData";
    }

    public String getSampleOutput(final SampleData sampleData)
    {
        String sampleDataStr = String.format("%d", sampleData.getSnvCount());

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

    private void writeSampleData(final SampleData sampleData)
    {
        try
        {
            final String sampleStr = String.format("%s,%d", sampleData.Id, sampleData.getSnvCount());

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

    private void loadRefSigContribPercentiles()
    {
        if(mConfig.RefSigContribData.isEmpty())
            return;

        populateRefSigContributions(mConfig.RefSigContribData, mRefCancerSigContribPercentiles);

        /*
        // cosmic_sig_subset.csv
        final GenericDataCollection sigsCollection = GenericDataLoader.loadFile(mConfig.RefSigContribData);
        mSnvSigNames.addAll(sigsCollection.getFieldNames());
        mSnvSignatures = DataUtils.createMatrixFromListData(sigsCollection.getData());
        mLeastSquaresFit = new LeastSquaresFit(mSnvSignatures.Rows, mSnvSignatures.Cols);
        */
    }

    public void close() { closeBufferedWriter(mSampleDataWriter); }

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

    /*
    private void fitSnvSignatures(final SampleData sampleData, int sampleCountsIndex)
    {
        if(mSnvSignatures == null)
            return;

        final double[] sampleCounts = mSampleCounts.getCol(sampleCountsIndex);

        mLeastSquaresFit.initialise(mSnvSignatures.getData(), sampleCounts);
        mLeastSquaresFit.solve();

        final double[] sigAllocs = mLeastSquaresFit.getContribs();

        double sampleTotal = sampleData.getSnvCount();

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
    */



}
