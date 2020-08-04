package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVERS_MIN_PREVALENCE_PERCENT;

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
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;

public class DriverAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final List<SampleDriverData> mSampleDrivers;
    private final Map<String,List<DriverPrevalence>> mCancerDriverPrevalence;
    private final SampleDataCache mSampleDataCache;
    private boolean mValidData;

    private final Set<String> mGeneList;
    private final Map<String,double[]> mGenePrevalenceTotals;

    private static final int POS_PREV = 0;
    private static final int NEG_PREV = 1;

    private BufferedWriter mSampleDataWriter;

    public DriverAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDrivers = Lists.newArrayList();
        mCancerDriverPrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mGeneList = Sets.newHashSet();
        mSampleDataCache = sampleDataCache;
        mSampleDataWriter = null;
        mValidData = false;

        if(config.DriverPrevFile.isEmpty() || config.SampleDriversFile.isEmpty())
            return;

        mValidData = true;
        loadPrevalenceData(config.DriverPrevFile);
        formGenePrevalenceTotals();
        loadSampleDrivers(config.SampleDriversFile);
        initialiseOutputFile();
    }

    public void processCohort()
    {
        if(!mValidData)
            return;

        mSampleDataCache.SampleDataList.forEach(x -> processSample(x));
        closeBufferedWriter(mSampleDataWriter);
    }

    public void processSample(final SampleData sample)
    {
        if(!mValidData)
            return;

        calculateSampleProbabilities(sample);
        writeSampleData(sample);
    }

    private void calculateSampleProbabilities(final SampleData sample)
    {
        final List<SampleDriverData> sampleDrivers = sample.getDrivers();

        if(sampleDrivers.isEmpty())
            return; // could still calculate a probability based on absence

        final Map<String,Double> cancerProbTotals = Maps.newHashMap();
        double allCancerProbTotal = 0;

        for(Map.Entry<String,List<DriverPrevalence>> entry : mCancerDriverPrevalence.entrySet())
        {
            final String cancerType = entry.getKey();

            if(cancerType.equals(CANCER_TYPE_UNKNOWN))
                continue;

            final List<DriverPrevalence> driverPrevalences = entry.getValue();

            double probabilityTotal = 1;

            Set<String> processedGenes = Sets.newHashSet();

            for(final SampleDriverData driverData : sampleDrivers)
            {
                if(processedGenes.contains(driverData.Gene))
                    continue;

                processedGenes.add(driverData.Gene);

                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverData.Gene);

                final DriverPrevalence driverPrev = driverPrevalences.stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(driverData.Gene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : DRIVERS_MIN_PREVALENCE_PERCENT;

                probabilityTotal *= driverPrevValue / genePrevTotals[POS_PREV];
            }

            allCancerProbTotal += probabilityTotal;
            cancerProbTotals.put(cancerType, probabilityTotal);
        }

        final Map<String,Double> driverCancerProbs = sample.getDriverCancerTypeProbs();

        for(Map.Entry<String,Double> entry : cancerProbTotals.entrySet())
        {
            double probability = entry.getValue() / allCancerProbTotal;

            if(probability > 0.05)
                driverCancerProbs.put(entry.getKey(), probability);
        }
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
                final DriverPrevalence data = DriverPrevalence.from(line);

                if(data == null)
                    continue;

                List<DriverPrevalence> dataList = mCancerDriverPrevalence.get(data.CancerType);
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
                final SampleDriverData data = SampleDriverData.from(line);

                if(data != null)
                {
                    SampleData sample = mSampleDataCache.SampleDataList.stream().filter(x -> x.Id.equals(data.SampleId)).findFirst().orElse(null);

                    if(sample != null)
                    {
                        sample.getDrivers().add(data);
                        mSampleDrivers.add(data);
                    }
                }
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
            for(Map.Entry<String,List<DriverPrevalence>> entry : mCancerDriverPrevalence.entrySet())
            {
                final DriverPrevalence driverPrevalence = entry.getValue().stream()
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

    public String getSampleOutput(final SampleData sample)
    {
        return String.format("NONE,0");
    }

    private void initialiseOutputFile()
    {
        try
        {
            mSampleDataWriter = createBufferedWriter(mConfig.formOutputFilename("DRIVERS"), false);

            mSampleDataWriter.write("SampleId,MatchedCancerType,Probability,DriverGenes");
            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample driver probability output: {}", e.toString());
        }
    }

    public void close() { closeBufferedWriter(mSampleDataWriter); }

    private void writeSampleData(final SampleData sample)
    {
        if(sample.getDrivers().isEmpty())
            return;

        try
        {
            String driverGenes = "";
            for(SampleDriverData driver : sample.getDrivers())
            {
                driverGenes = appendStr(driverGenes, driver.Gene, ';');
            }

            if(sample.getDriverCancerTypeProbs().isEmpty())
            {
                mSampleDataWriter.write(String.format("%s,NONE,0,%s",
                        sample.Id, driverGenes));
                mSampleDataWriter.newLine();
                return;
            }

            for(Map.Entry<String,Double> driverCancerProbs : sample.getDriverCancerTypeProbs().entrySet())
            {
                mSampleDataWriter.write(String.format("%s,%s,%.4f,%s",
                        sample.Id, driverCancerProbs.getKey(), driverCancerProbs.getValue(), driverGenes));
                mSampleDataWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample driver probability output: {}", e.toString());
        }
    }

/*

            final String cancerType = entry.getKey();

            if(cancerType.equals(CANCER_TYPE_UNKNOWN))
                continue;

            final List<DriverPrevalence> driverPrevalences = entry.getValue();

            double probabilityTotal = 1;

            for(String gene : mGeneList)
            {
                final double[] genePrevTotals = mGenePrevalenceTotals.get(gene);

                final DriverPrevalence driverPrev = driverPrevalences.stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : DRIVERS_MIN_PREVALENCE_PERCENT;

                if(hasGeneMatch(gene, sampleDrivers))
                {
                    probabilityTotal *= driverPrevValue / genePrevTotals[POS_PREV];
                }
                else
                {
                    probabilityTotal *= (1 - driverPrevValue) / genePrevTotals[NEG_PREV];
                }
            }

            final List<String> matchedGenes = Lists.newArrayList();

            for(final DriverPrevalence driverPrev : driverPrevalences)
            {
                // for now skip mutation-type prevalence
                if(!driverPrev.isTypeAll())
                    continue;

                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverPrev.Gene);

                if(hasGeneMatch(driverPrev.Gene, sampleDrivers))
                {
                    matchedGenes.add(driverPrev.Gene);
                    probabilityTotal *= driverPrev.Prevalence / genePrevTotals[POS_PREV];
                }
                else
                {
                    probabilityTotal *= (1 - driverPrev.Prevalence) / genePrevTotals[NEG_PREV];
                }
            }

            for(final SampleDriverData driverData : sampleDrivers)
            {
                if(matchedGenes.contains(driverData.Gene))
                    continue;

                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverData.Gene);

                probabilityTotal *= DRIVERS_MIN_PREVALENCE_PERCENT / genePrevTotals[POS_PREV];
            }

            allCancerProbTotal += probabilityTotal;
            cancerProbTotals.put(cancerType, probabilityTotal);
 */
}
