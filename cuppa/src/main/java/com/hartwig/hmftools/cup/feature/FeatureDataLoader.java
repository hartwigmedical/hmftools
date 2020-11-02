package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSS;
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
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.cup.sigs.SomaticDataLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantFusionDAO;

public class FeatureDataLoader
{
    public static boolean loadDriversFromCohortFile(final String filename, final Map<String,List<SampleFeatureData>> sampleDrivers)
    {
        if(filename == null)
            return true;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

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

    public static boolean loadFeaturesFromFile(
            final String sampleId, final String sampleDataDir, final String sampleVcfFile,
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        try
        {
            String viralInsertFilename = LinxViralInsertion.generateFilename(sampleDataDir, sampleId);

            final List<LinxViralInsertion> viralInserts = Files.exists(Paths.get(viralInsertFilename)) ?
                    LinxViralInsertion.read(viralInsertFilename) : Lists.newArrayList();

            final String fusionsFilename = LinxFusion.generateFilename(sampleDataDir, sampleId);

            final List<LinxFusion> fusions = Files.exists(Paths.get(fusionsFilename)) ?
                    LinxFusion.read(fusionsFilename) : Lists.newArrayList();

            final String driverCatalogFilename = DriverCatalogFile.generateFilename(sampleDataDir, sampleId);

            final List<DriverCatalog> drivers = Files.exists(Paths.get(driverCatalogFilename)) ?
                    DriverCatalogFile.read(driverCatalogFilename) : Lists.newArrayList();

            final List<SomaticVariant> indels = checkIndels(sampleId, sampleDataDir, null) ?
                    SomaticDataLoader.loadSomaticVariants(sampleId, sampleVcfFile, Lists.newArrayList(VariantType.INDEL)) : null;

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, viralInserts, indels);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load drivers, fusion and virus data files: {}", e.toString());
            return false;
        }

        return true;
    }

    private static boolean checkIndels(final String sampleId, final String sampleDataDir, final DatabaseAccess dbAccess)
    {
        PurityContext purityContext = null;

        if(dbAccess != null)
        {
            purityContext = dbAccess.readPurityContext(sampleId);
        }
        else
        {
            try
            {
                purityContext = PurityContextFile.read(sampleDataDir, sampleId);
            }
            catch (Exception e)
            {
                CUP_LOGGER.error("sample({}) failed to load purity file( from dir{}): {}",
                        sampleId, sampleDataDir, e.toString());
                return false;
            }
        }

        return purityContext.microsatelliteStatus() == MSS;
    }

    public static boolean loadFeaturesFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        if(dbAccess == null)
            return false;

        int i = 0;
        int nextLog = 100;

        for(final String sampleId : sampleIds)
        {
            final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);
            final List<LinxFusion> fusions = dbAccess.readFusions(sampleId);
            final List<LinxViralInsertion> viralInserts = dbAccess.readViralInsertions(sampleId);

            final List<SomaticVariant> indels = checkIndels(sampleId, null, dbAccess) ?
                    dbAccess.readSomaticVariants(sampleId, VariantType.INDEL) : null;

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, viralInserts, indels);

            ++i;
            if(i >= nextLog)
            {
                nextLog += 100;
                CUP_LOGGER.info("loaded {} sample feature data sets", i);
            }
        }

        return true;
    }

    private static void mapFeatureData(
            final String sampleId, final Map<String,List<SampleFeatureData>> sampleDrivers, final List<DriverCatalog> drivers,
            final List<LinxFusion> fusions, final List<LinxViralInsertion> viralInserts, final List<SomaticVariant> indels)
    {
        final List<SampleFeatureData> featuresList = Lists.newArrayList();

        if(drivers != null)
        {
            for(final DriverCatalog driver : drivers)
            {
                SampleFeatureData feature = new SampleFeatureData(sampleId, driver.gene(), DRIVER, driver.driverLikelihood());

                if(driver.driver() == DriverType.AMP)
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_AMP);
                else if(driver.driver() == DriverType.DEL)
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_DEL);

                feature.ExtraInfo.put(DRIVER_CHROMOSOME, driver.chromosome());

                featuresList.add(feature);
            }
        }

        if(fusions != null)
        {
            for(final LinxFusion fusion : fusions)
            {
                if(!fusion.reported())
                    continue;

                if(fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString())
                || fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_3.toString()))
                {
                    if(fusion.likelihood() != FusionLikelihoodType.HIGH)
                        continue;

                    final String[] genes = fusion.name().split("_");

                    final String fusionName = fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString()) ?
                            genes[0] + "_PROM5" : genes[1] + "_PROM3";

                    featuresList.add(new SampleFeatureData(sampleId, fusionName, FeatureType.FUSION, 1));
                }
                else
                {
                    featuresList.add(new SampleFeatureData(sampleId, fusion.name(), FeatureType.FUSION, 1));
                }
            }
        }

        if(viralInserts != null)
        {
            final List<SampleFeatureData> viralInsertDataList = viralInserts.stream()
                    .map(x -> new SampleFeatureData(sampleId, fromVirusName(x.VirusName).toString(), FeatureType.VIRUS, 1))
                    .filter(x -> !x.Name.equals(OTHER.toString()))
                    .collect(Collectors.toList());

            featuresList.addAll(viralInsertDataList);
        }

        if(indels != null)
        {
            final List<SampleFeatureData> indelFeatures = indels.stream()
                    .filter(x -> x.gene().equals(INDEL_ALB) || x.gene().equals(INDEL_SFTPB) || x.gene().equals(INDEL_SLC34A2))
                    .map(x -> new SampleFeatureData(sampleId, String.format("INDEL_%s", x.gene()), FeatureType.INDEL, 1))
                    .filter(x -> !x.Name.equals(OTHER.toString()))
                    .collect(Collectors.toList());

            featuresList.addAll(indelFeatures);
        }

        sampleDrivers.put(sampleId, featuresList);
    }

    private static final String INDEL_ALB = "ALB";
    private static final String INDEL_SFTPB = "SFTPB";
    private static final String INDEL_SLC34A2 = "SLC34A2";

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
