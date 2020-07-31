package com.hartwig.hmftools.sig_analyser.cup;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.CUP_LOGGER;
import static com.hartwig.hmftools.sig_analyser.cup.CupConstants.DRIVERS_MIN_PREVALENCE_PERCENT;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class CupDrivers
{
    private final CupConfig mConfig;
    private final List<CupSampleDriverData> mSampleDrivers;
    private final Map<String,List<CupDriverPrevalence>> mCancerDriverPrevalence;
    private final List<CupSampleData> mSampleDataList;

    private final Set<String> mGeneList;
    private final Map<String,double[]> mGenePrevalenceTotals;

    private static final int POS_PREV = 0;
    private static final int NEG_PREV = 1;

    private BufferedWriter mSampleDataWriter;

    public CupDrivers(final CupConfig config, final List<CupSampleData> sampleDataList)
    {
        mConfig = config;
        mSampleDrivers = Lists.newArrayList();
        mCancerDriverPrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mGeneList = Sets.newHashSet();
        mSampleDataList = sampleDataList;
        mSampleDataWriter = null;

        loadPrevalenceData(config.DriverPrevFile);
        formGenePrevalenceTotals();
        loadSampleDrivers(config.SampleDriversFile);
    }

    public void run()
    {
        mSampleDataList.forEach(x -> calculateSampleProbabilities(x));
    }

    private void calculateSampleProbabilities(final CupSampleData sampleData)
    {
        final List<CupSampleDriverData> sampleDrivers = sampleData.getDrivers();

        // DRIVERS_MIN_PREVALENCE_PERCENT

        final Map<String,Double> cancerProbTotals = Maps.newHashMap();
        double allCancerProbTotal = 0;

        for(Map.Entry<String,List<CupDriverPrevalence>> entry : mCancerDriverPrevalence.entrySet())
        {
            final String cancerType = entry.getKey();

            final List<CupDriverPrevalence> driverPrevalences = entry.getValue();

            double probabilityTotal = 0;

            final List<String> matchedGenes = Lists.newArrayList();

            for(final CupDriverPrevalence driverPrev : driverPrevalences)
            {
                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverPrev.Gene);

                if(hasGeneMatch(driverPrev.Gene, sampleDrivers, false))
                {
                    matchedGenes.add(driverPrev.Gene);
                    probabilityTotal += driverPrev.Prevalence / genePrevTotals[POS_PREV];
                }
                else
                {
                    probabilityTotal += (1 - driverPrev.Prevalence) / genePrevTotals[NEG_PREV];
                }
            }

            for(final CupSampleDriverData driverData : sampleDrivers)
            {
                if(matchedGenes.contains(driverData.Gene))
                    continue;

                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverData.Gene);

                probabilityTotal += DRIVERS_MIN_PREVALENCE_PERCENT / genePrevTotals[POS_PREV];
            }

            allCancerProbTotal += probabilityTotal;
            cancerProbTotals.put(cancerType, probabilityTotal);
        }

        final Map<String,Double> driverCancerProbs = sampleData.getDriverCancerTypeProbs();

        for(Map.Entry<String,Double> entry : cancerProbTotals.entrySet())
        {
            driverCancerProbs.put(entry.getKey(), entry.getValue() / allCancerProbTotal);
        }
    }

    private boolean hasGeneMatch(final String gene, final List<CupSampleDriverData> sampleDrivers, boolean useSpecificTypes)
    {
        return sampleDrivers.stream().filter(x -> useSpecificTypes != x.isTypeAll()).anyMatch(x -> x.Gene.equals(gene));
    }

    private void loadPrevalenceData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final CupDriverPrevalence data = CupDriverPrevalence.from(line);

                if(data == null)
                    continue;

                List<CupDriverPrevalence> dataList = mCancerDriverPrevalence.get(data.CancerType);
                if(dataList == null)
                {
                    mCancerDriverPrevalence.put(data.CancerType, Lists.newArrayList(data));
                }
                else
                {
                    dataList.add(data);
                }

                mGeneList.add(data.Gene);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read driver prevalence data file({}): {}", filename, e.toString());
        }
    }

    private void loadSampleDrivers(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final CupSampleDriverData data = CupSampleDriverData.from(line);

                if(data != null)
                    mSampleDrivers.add(data);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample driver data file({}): {}", filename, e.toString());
        }
    }

    private void formGenePrevalenceTotals()
    {
        for(final String gene : mGeneList)
        {
            double[] genePrevTotals = new double[NEG_PREV+1];
            mGenePrevalenceTotals.put(gene, genePrevTotals);
            for(Map.Entry<String,List<CupDriverPrevalence>> entry : mCancerDriverPrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final CupDriverPrevalence driverPrevalence = entry.getValue().stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(driverPrevalence != null)
                {
                    genePrevTotals[POS_PREV] += driverPrevalence.Prevalence;
                    genePrevTotals[NEG_PREV] += 1 - driverPrevalence.Prevalence;
                }
                else
                {
                    genePrevTotals[POS_PREV] += DRIVERS_MIN_PREVALENCE_PERCENT;
                    genePrevTotals[NEG_PREV] += 1 - DRIVERS_MIN_PREVALENCE_PERCENT;
                }
            }
        }
    }

    public String getHeader()
    {
        return "TopDriverCancerType,TopDriverProb";
    }

    public String getSampleOutput(final CupSampleData sampleData)
    {
        return String.format(",,");
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
            /*
            final String sampleStr = String.format("%s,%d",
                    sampleData.SampleId, mSampleTotals[sampleData.index()]);


            for(Map.Entry<String,Double> cancerCssData : cssDataMap.entrySet())
            {
                mSampleDataWriter.write(String.format("%s,%.4f,%s,%.4f,%s",
                        sampleStr, totalWeightedCss, cancerCssData.getKey(), cancerCssData.getValue(), sigDataStr));
                mSampleDataWriter.newLine();
            }

            */

            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SNV sample CSS output: {}", e.toString());
        }
    }


}
