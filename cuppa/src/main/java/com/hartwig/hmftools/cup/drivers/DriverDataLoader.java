package com.hartwig.hmftools.cup.drivers;

import static java.lang.Math.max;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_LIKELIHOOD_THRESHOLD;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_MIN_PREVALENCE;
import static com.hartwig.hmftools.cup.common.CupConstants.FUSION_MIN_PREVALENCE;
import static com.hartwig.hmftools.cup.drivers.DriverType.DRIVER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverDataLoader
{
    public static void loadDriversFromCohortFile(final String filename, final Map<String,List<SampleDriverData>> sampleDrivers)
    {
        if(filename == null)
            return;

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
                    List<SampleDriverData> drivers = sampleDrivers.get(driverData.SampleId);

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
        } catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample driver data file({}): {}", filename, e.toString());
        }
    }

    public static void loadDriversFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,List<SampleDriverData>> sampleDrivers)
    {
        if(dbAccess == null)
            return;

        for(final String sampleId : sampleIds)
        {
            final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);
            if(drivers != null)
            {
                final List<SampleDriverData> driverDataList = drivers.stream()
                        .filter(x -> x.driverLikelihood() >= DRIVER_LIKELIHOOD_THRESHOLD)
                        .map(x -> new SampleDriverData(sampleId, x.gene(), DRIVER))
                        .collect(Collectors.toList());

                sampleDrivers.put(sampleId, driverDataList);
            }

            final List<LinxFusion> fusions =dbAccess.readFusions(sampleId);

            if(fusions != null)
            {
                final List<SampleDriverData> fusionDataList = fusions.stream()
                        .filter(x -> x.reported())
                        .map(x -> new SampleDriverData(sampleId, x.name(), DriverType.FUSION))
                        .collect(Collectors.toList());

                sampleDrivers.put(sampleId, fusionDataList);
            }
        }
    }

    public static void loadRefPrevalenceData(
            final String filename, final Map<String,DriverPrevCounts> genePrevalenceTotals,
            final Map<String,List<DriverPrevData>> cancerDriverPrevalence)
    {
        if(filename == null || filename.isEmpty())
            return;

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

                DriverPrevCounts genePrevTotals = genePrevalenceTotals.get(prevData.Gene);

                if(genePrevTotals == null)
                {
                    genePrevTotals = new DriverPrevCounts();
                    genePrevalenceTotals.put(prevData.Gene, genePrevTotals);
                    genePrevTotals.MinPrevalence = prevData.Type == DRIVER ? DRIVER_MIN_PREVALENCE : FUSION_MIN_PREVALENCE;
                }

                genePrevTotals.MaxPrevalence = max(genePrevTotals.MaxPrevalence, prevData.Prevalence);

                final List<DriverPrevData> dataList = cancerDriverPrevalence.get(prevData.CancerType);
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
        }
    }



}
