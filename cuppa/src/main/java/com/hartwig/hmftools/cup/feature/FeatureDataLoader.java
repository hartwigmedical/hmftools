package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_CHROMOSOME;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.cup.feature.ViralInsertionType.OTHER;
import static com.hartwig.hmftools.cup.feature.ViralInsertionType.fromVirusName;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class FeatureDataLoader
{
    public static boolean loadDriversFromCohortFile(final String filename, final Map<String,List<SampleFeatureData>> sampleDrivers)
    {
        if(filename == null)
            return true;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final SampleFeatureData driverData = SampleFeatureData.from(line);

                if(driverData != null)
                {
                    List<SampleFeatureData> drivers = sampleDrivers.get(driverData.SampleId);

                    if(drivers == null)
                    {
                        sampleDrivers.put(driverData.SampleId, Lists.newArrayList(driverData));
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
            return false;
        }

        return true;
    }

    public static boolean loadDriversFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,List<SampleFeatureData>> sampleDrivers)
    {
        if(dbAccess == null)
            return false;

        for(final String sampleId : sampleIds)
        {
            final List<SampleFeatureData> featuresList = Lists.newArrayList();

            final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);
            if(drivers != null)
            {

                for(final DriverCatalog driver : drivers)
                {
                    if(driver.driverLikelihood() == 0)
                        continue;

                    SampleFeatureData feature = new SampleFeatureData(sampleId, driver.gene(), DRIVER, driver.driverLikelihood());

                    if(driver.driver() == DriverType.AMP)
                        feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_AMP);
                    else if(driver.driver() == DriverType.DEL)
                        feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_DEL);

                    feature.ExtraInfo.put(DRIVER_CHROMOSOME, driver.chromosome());

                    featuresList.add(feature);
                }
            }

            final List<LinxFusion> fusions = dbAccess.readFusions(sampleId);

            if(fusions != null)
            {
                final List<SampleFeatureData> fusionDataList = fusions.stream()
                        .filter(x -> x.reported())
                        .map(x -> new SampleFeatureData(sampleId, x.name(), FeatureType.FUSION, 1))
                        .collect(Collectors.toList());

                featuresList.addAll(fusionDataList);
            }

            final List<LinxViralInsertion> viralInserts = dbAccess.readViralInsertions(sampleId);

            if(viralInserts != null)
            {
                final List<SampleFeatureData> viralInsertDataList = viralInserts.stream()
                        .map(x -> new SampleFeatureData(sampleId, fromVirusName(x.VirusName).toString(), FeatureType.VIRUS, 1))
                        .filter(x -> !x.Gene.equals(OTHER.toString()))
                        .collect(Collectors.toList());

                featuresList.addAll(viralInsertDataList);
            }

            sampleDrivers.put(sampleId, featuresList);
        }

        return true;
    }

    public static boolean loadRefPrevalenceData(
            final String filename, final Map<String,FeaturePrevCounts> genePrevalenceTotals,
            final Map<String,List<FeaturePrevData>> cancerDriverPrevalence)
    {
        if(filename == null || filename.isEmpty())
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                final FeaturePrevData prevData = FeaturePrevData.from(line);

                if(prevData == null)
                    continue;

                FeaturePrevCounts genePrevTotals = genePrevalenceTotals.get(prevData.Gene);

                if(genePrevTotals == null)
                {
                    genePrevTotals = new FeaturePrevCounts();
                    genePrevalenceTotals.put(prevData.Gene, genePrevTotals);
                    // genePrevTotals.MinPrevalence = prevData.Type == DRIVER ? DRIVER_MIN_PREVALENCE : FUSION_MIN_PREVALENCE;
                }

                genePrevTotals.MaxPrevalence = max(genePrevTotals.MaxPrevalence, prevData.Prevalence);

                final List<FeaturePrevData> dataList = cancerDriverPrevalence.get(prevData.CancerType);
                if(dataList == null)
                {
                    cancerDriverPrevalence.put(prevData.CancerType, Lists.newArrayList(prevData));
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
            return false;
        }

        return true;
    }

    public static boolean loadRefCancerFeatureAvg(final String filename, final Map<String,Double> cancerFeatureAvgs)
    {
        if(filename == null || filename.isEmpty())
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM);
                final String cancerType = items[0];
                double average = Double.parseDouble(items[1]);
                cancerFeatureAvgs.put(cancerType, average);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read feature averages data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

}
