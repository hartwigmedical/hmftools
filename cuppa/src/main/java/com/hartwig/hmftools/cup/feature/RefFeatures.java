package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_FEATURE_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_DRIVER_AVG;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_FEATURE_PREV;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefFeatureOverrides;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.parseFileSet;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.ref.RefClassifier;

public class RefFeatures implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<FeaturePrevData>> mFeatureOverrides;

    public RefFeatures(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mFeatureOverrides = Maps.newHashMap();
        loadRefFeatureOverrides(mConfig.FeatureOverrideFile, mFeatureOverrides);
    }

    public CategoryType categoryType() { return FEATURE; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.DbAccess != null || !config.CohortFeaturesFile.isEmpty() || !config.SampleFeaturesDir.isEmpty();
    }

    public void buildRefDataSets()
    {
        if(mConfig.CohortFeaturesFile.isEmpty() && mConfig.DbAccess == null && mConfig.SampleFeaturesDir.isEmpty())
            return;

        CUP_LOGGER.info("building feature reference data");

        final Map<String,List<SampleFeatureData>> sampleFeaturesMap = Maps.newHashMap();
        loadSampleData(sampleFeaturesMap);

        writeCohortData(sampleFeaturesMap);

        final Map<String,Map<String,Double>> cancerFeatureCounts = Maps.newHashMap();
        final Map<String,FeatureType> featureTypes = Maps.newHashMap();

        final Map<String,List<Double>> driversPerSampleMap = Maps.newHashMap();
        final List<Double> panCancerDriversPerSample = Lists.newArrayList();

        assignFeatures(sampleFeaturesMap, cancerFeatureCounts, featureTypes, driversPerSampleMap, panCancerDriversPerSample);

        writeFeaturePrevalenceFile(cancerFeatureCounts, featureTypes);

        writeAverageDriversFile(driversPerSampleMap, panCancerDriversPerSample);
    }

    private void loadSampleData(final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        if(!mConfig.CohortFeaturesFile.isEmpty())
        {
            final Map<String,List<SampleFeatureData>> allSampleFeatures = Maps.newHashMap();

            final List<String> files = parseFileSet(mConfig.CohortFeaturesFile);
            files.forEach(x -> loadFeaturesFromCohortFile(x, allSampleFeatures));

            // extract only reference sample data
            for(Map.Entry<String,List<SampleFeatureData>> entry : allSampleFeatures.entrySet())
            {
                String sampleId = entry.getKey();

                if(mSampleDataCache.hasRefSample(sampleId))
                {
                    sampleFeaturesMap.put(sampleId, entry.getValue());
                }
            }
        }
        else if(mConfig.DbAccess != null)
        {
            loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.refSampleIds(false), sampleFeaturesMap);
        }
        else
        {
            // load from per-sample files
            for(int i = 0; i < mSampleDataCache.RefSampleDataList.size(); ++i)
            {
                SampleData sample = mSampleDataCache.RefSampleDataList.get(i);

                final String sampleDataDir = formSamplePath(mConfig.SampleFeaturesDir, sample.Id);
                final String somaticVcf = formSamplePath(mConfig.SampleSomaticVcf, sample.Id);

                if(!loadFeaturesFromFile(sample.Id, sampleDataDir, somaticVcf, sampleFeaturesMap))
                    break;

                if(i > 0 && (i % 100) == 0)
                {
                    CUP_LOGGER.debug("processed {} sample feature files", i);
                }
            }
        }
    }

    private void assignFeatures(
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap,
            final Map<String,Map<String,Double>> cancerFeatureCounts, final Map<String,FeatureType> featureTypes,
            final Map<String,List<Double>> driversPerSampleMap, final List<Double> panCancerDriversPerSample)
    {
        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            final String sampleId = sampleEntry.getKey();
            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

            if(!isKnownCancerType(cancerType))
                continue;

            final Set<String> features = Sets.newHashSet();
            double driverTotal = 0;

            for(SampleFeatureData feature : sampleEntry.getValue())
            {
                if(feature.Likelihood == 0)
                    continue;

                if(feature.Type == DRIVER)
                    driverTotal += feature.Likelihood;

                if(features.contains(feature.Name))
                    continue;

                features.add(feature.Name);
                featureTypes.put(feature.Name, feature.Type);

                double likelihoodTotal = sampleEntry.getValue().stream().filter(x -> x.Name.equals(feature.Name)).mapToDouble(x -> x.Likelihood).sum();
                likelihoodTotal = min(likelihoodTotal, 1);

                Map<String,Double> featureCounts = cancerFeatureCounts.get(cancerType);
                if(featureCounts == null)
                {
                    featureCounts = Maps.newHashMap();
                    cancerFeatureCounts.put(cancerType, featureCounts);
                }

                Double featureCount = featureCounts.get(feature.Name);

                if(featureCount == null)
                    featureCounts.put(feature.Name, likelihoodTotal);
                else
                    featureCounts.put(feature.Name, likelihoodTotal + featureCount);
            }

            List<Double> sampleDrivers = driversPerSampleMap.get(cancerType);
            if(sampleDrivers == null)
            {
                sampleDrivers = Lists.newArrayList(driverTotal);
                driversPerSampleMap.put(cancerType, sampleDrivers);
            }
            else
            {
                sampleDrivers.add(driverTotal);
            }

            panCancerDriversPerSample.add(driverTotal);
        }
    }

    private void writeFeaturePrevalenceFile(final Map<String,Map<String,Double>> cancerFeatureCounts, final Map<String,FeatureType> featureTypes)
    {
        // output: CancerType,Gene,Type,SamplePerc
        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_FEATURE_PREV;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("CancerType,Feature,Type,SamplePerc");
            writer.newLine();

            for(Map.Entry<String,Map<String,Double>> cancerEntry : cancerFeatureCounts.entrySet())
            {
                final String cancerType = cancerEntry.getKey();
                int cancerSamples = mSampleDataCache.getCancerSampleCount(cancerType);

                for(Map.Entry<String,Double> featureEntry : cancerEntry.getValue().entrySet())
                {
                    final String feature = featureEntry.getKey();

                    // check for any feature-cancer overrides
                    final List<FeaturePrevData> overrides = mFeatureOverrides.get(feature);

                    if(overrides == null)
                    {
                        double sampleTotal = featureEntry.getValue();
                        double prevalence = sampleTotal / cancerSamples;
                        FeatureType featureType = featureTypes.get(feature);

                        writer.write(String.format("%s,%s,%s,%.6f", cancerType, feature, featureType, prevalence));
                        writer.newLine();
                    }
                    else
                    {
                        // remove this entry since unobserved features will then be written at the end
                        final FeaturePrevData override = overrides.stream().filter(x -> x.CancerType.equals(cancerType)).findFirst().orElse(null);

                        if(override != null)
                        {
                            writer.write(String.format("%s,%s,%s,%.6f", cancerType, feature, override.Type, override.Prevalence));
                            writer.newLine();

                            overrides.remove(override);
                        }
                    }
                }
            }

            for(List<FeaturePrevData> overrides : mFeatureOverrides.values())
            {
                if(overrides.isEmpty())
                    continue;

                for(FeaturePrevData override : overrides)
                {
                    writer.write(String.format("%s,%s,%s,%.6f", override.CancerType, override.Name, override.Type, override.Prevalence));
                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref feature prevalence data output: {}", e.toString());
        }
    }

    private void writeAverageDriversFile(final Map<String,List<Double>> driversPerSampleMap, final List<Double> panCancerDriversPerSample)
    {
        // CancerType,AvgFeatures
        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_DRIVER_AVG;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("CancerType,AvgDrivers");
            writer.newLine();

            for(Map.Entry<String,List<Double>> entry : driversPerSampleMap.entrySet())
            {
                final String cancerType = entry.getKey();
                final List<Double> drivers = entry.getValue();

                double driverTotal = drivers.stream().mapToDouble(x -> x).sum();
                int cancerSamples = mSampleDataCache.getCancerSampleCount(cancerType);
                double avgDrivers = driverTotal / cancerSamples;

                CUP_LOGGER.debug("cancerType({}) samples({}) driverTotal({}) avgPerSample({})",
                        cancerType, cancerSamples, String.format("%.2f", driverTotal), String.format("%.4f", avgDrivers));

                writer.write(String.format("%s,%.2f", cancerType, avgDrivers));
                writer.newLine();
            }

            double driverTotal = panCancerDriversPerSample.stream().mapToDouble(x -> x).sum();

            long knownRefSampleCount = mSampleDataCache.RefSampleDataList.stream().filter(x -> isKnownCancerType(x.cancerType())).count();
            double avgDrivers = driverTotal / knownRefSampleCount;

            CUP_LOGGER.debug("pan-cancer driverTotal({}) avgPerSample({}) knownRefSampleCount({})",
                    String.format("%.2f", driverTotal), String.format("%.4f", avgDrivers), knownRefSampleCount);

            writer.write(String.format("ALL,%.2f", avgDrivers));
            writer.newLine();

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref driver averages file: {}", e.toString());
        }
    }

    private void writeCohortData(final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        if(!mConfig.WriteCohortFiles)
            return;

        final String filename = mConfig.OutputDir + COHORT_REF_FILE_FEATURE_DATA;
        if(Files.exists(Paths.get(filename)))
        {
            CUP_LOGGER.warn("not over-writing cohort feature reference file({})", filename);
            return;
        }

        CUP_LOGGER.info("writing cohort feature reference data");

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(SampleFeatureData.header());
            writer.newLine();

            for(Map.Entry<String,List<SampleFeatureData>> entry : sampleFeaturesMap.entrySet())
            {
                final String sampleId = entry.getKey();

                for(SampleFeatureData feature : entry.getValue())
                {
                    writer.write(String.format("%s,%s", sampleId, feature.toCsv()));
                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write feature cohort data output: {}", e.toString());
        }
    }
}
