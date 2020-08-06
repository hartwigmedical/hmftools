package com.hartwig.hmftools.cup.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.DRIVER;
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
import com.hartwig.hmftools.cup.common.SampleResult;

public class DriverAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final Map<String,List<SampleDriverData>> mSampleDrivers;
    private final Map<String,List<DriverPrevData>> mCancerDriverPrevalence;
    private final SampleDataCache mSampleDataCache;
    private boolean mValidData;

    private final Map<String,double[]> mGenePrevalenceTotals;

    private static final int POS_PREV = 0;
    private static final int NEG_PREV = 1;
    private static final int MAX_PREV = 2;

    private BufferedWriter mSampleDataWriter;

    public DriverAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDrivers = Maps.newHashMap();
        mCancerDriverPrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mSampleDataWriter = null;
        mValidData = false;

        if(config.RefDriverPrevFile.isEmpty() || config.SampleDriversFile.isEmpty())
            return;

        mValidData = true;
        loadPrevalenceData(config.RefDriverPrevFile);
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

    public final List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(!mValidData)
            return results;

        final List<SampleDriverData> sampleDrivers = mSampleDrivers.get(sample.Id);

        if(sampleDrivers == null || sampleDrivers.isEmpty())
            return results;

        addDriverPrevalence(sample, sampleDrivers, results);
        calcCancerTypeProbability(sample, sampleDrivers, results);

        // writeSampleData(sample);

        return results;
    }

    private void addDriverPrevalence(final SampleData sample, final List<SampleDriverData> sampleDrivers, final List<SampleResult> results)
    {
        for(final SampleDriverData driverData : sampleDrivers)
        {
            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,List<DriverPrevData>> entry : mCancerDriverPrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final List<DriverPrevData> driverPrevalences = entry.getValue();

                final DriverPrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(driverData.Gene)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.Prevalence : 0);
            }

            SampleResult result = new SampleResult(sample.Id, DRIVER, driverData.Type, driverData.Gene, cancerTypeValues);
            results.add(result);
        }
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<SampleDriverData> sampleDrivers, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String, Double> cancerProbTotals = Maps.newHashMap();
        double allCancerProbTotal = 0;

        final Set<String> genes = Sets.newHashSet();
        sampleDrivers.forEach(x -> genes.add(x.Gene));

        for(Map.Entry<String, List<DriverPrevData>> entry : mCancerDriverPrevalence.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<DriverPrevData> driverPrevalences = entry.getValue();

            double probabilityTotal = 1;

            for(final String driverGene : genes)
            {
                final double[] genePrevTotals = mGenePrevalenceTotals.get(driverGene);
                double minPrevalence = min(MIN_PREV_PERC_OF_MAX * genePrevTotals[MAX_PREV], MIN_PREVALENCE);

                final DriverPrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(driverGene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : minPrevalence;
                probabilityTotal *= driverPrevValue / genePrevTotals[POS_PREV];
            }

            cancerProbTotals.put(cancerType, probabilityTotal);
            allCancerProbTotal += probabilityTotal;
        }

        final Map<String,Double> cancerTypeValues = Maps.newHashMap();

        for(Map.Entry<String,Double> entry : cancerProbTotals.entrySet())
        {
            double probability = entry.getValue() / allCancerProbTotal;
            cancerTypeValues.put(entry.getKey(), probability);
        }

        SampleResult result = new SampleResult(sample.Id, CLASSIFIER, "DRIVER", "COMBINED", cancerTypeValues);
        results.add(result);
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
                final DriverPrevData prevData = DriverPrevData.from(line);

                if(prevData == null)
                    continue;

                double[] genePrevTotals = mGenePrevalenceTotals.get(prevData.Gene);

                if(genePrevTotals == null)
                {
                    genePrevTotals = new double[MAX_PREV + 1];
                    mGenePrevalenceTotals.put(prevData.Gene, genePrevTotals);
                }

                genePrevTotals[MAX_PREV] = max(genePrevTotals[MAX_PREV], prevData.Prevalence);

                final List<DriverPrevData> dataList = mCancerDriverPrevalence.get(prevData.CancerType);
                if(dataList == null)
                {
                    mCancerDriverPrevalence.put(prevData.CancerType, Lists.newArrayList(prevData));
                }
                else
                {
                    dataList.add(prevData);
                }
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
                final SampleDriverData driverData = SampleDriverData.from(line);

                if(driverData != null)
                {
                    List<SampleDriverData> drivers = mSampleDrivers.get(driverData.SampleId);

                    if(drivers == null)
                    {
                        mSampleDrivers.put(driverData.SampleId, Lists.newArrayList(driverData));
                    }
                    else
                    {
                        drivers.add(driverData);
                    }
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample driver data file({}): {}", filename, e.toString());
        }
    }

    private static final double MIN_PREV_PERC_OF_MAX = 0.1;
    private static final double MIN_PREVALENCE = 0.01;

    private void formGenePrevalenceTotals()
    {
        for(Map.Entry<String,double[]> geneEntry : mGenePrevalenceTotals.entrySet())
        {
            final String gene = geneEntry.getKey();
            final double[] genePrevTotals = geneEntry.getValue();

            double minPrevalence = min(MIN_PREV_PERC_OF_MAX * genePrevTotals[MAX_PREV], MIN_PREVALENCE);

            for(Map.Entry<String,List<DriverPrevData>> cancerEntry : mCancerDriverPrevalence.entrySet())
            {
                final DriverPrevData driverPrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(driverPrevData != null)
                {
                    genePrevTotals[POS_PREV] += driverPrevData.Prevalence;
                    genePrevTotals[NEG_PREV] += 1 - driverPrevData.Prevalence;
                }
                else
                {
                    genePrevTotals[POS_PREV] += minPrevalence;
                    genePrevTotals[NEG_PREV] += 1 - minPrevalence;
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
