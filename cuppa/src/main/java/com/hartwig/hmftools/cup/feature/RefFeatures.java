package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FEATURE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_DRIVER_AVG;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_FEATURE_PREV;
import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CupConstants.loadKnownMutations;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefFeatureOverrides;
import static com.hartwig.hmftools.cup.feature.FeaturePrevData.TYPE_NAME_DELIM;
import static com.hartwig.hmftools.cup.feature.FeaturePrevData.featureTypeName;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.convertDriverAmps;
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
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.ref.RefClassifier;

public class RefFeatures implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final List<FeaturePrevData> mFeatureOverrides;

    public RefFeatures(final RefDataConfig config, final SampleDataCache sampleDataCache, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        loadKnownMutations(mConfig.RefGenVersion);

        mSampleDataCache = sampleDataCache;

        mFeatureOverrides = loadRefFeatureOverrides(mConfig.FeatureOverrideFile);
    }

    public CategoryType categoryType() { return FEATURE; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        if(config.Categories.contains(FEATURE))
            return true;

        return config.DbAccess != null || !config.CohortFeaturesFile.isEmpty() || !config.LinxDir.isEmpty();
    }

    @Override
    public boolean buildRefDataSets()
    {
        if(mConfig.CohortFeaturesFile.isEmpty() && mConfig.DbAccess == null && mConfig.LinxDir.isEmpty())
        {
            CUP_LOGGER.error("feature ref builder missing DB config or directories");
            return false;
        }

        CUP_LOGGER.info("building feature reference data");

        final Map<String,List<SampleFeatureData>> sampleFeaturesMap = Maps.newHashMap();

        if(!loadSampleData(sampleFeaturesMap))
            return false;

        writeCohortData(sampleFeaturesMap);

        convertDriverAmps(sampleFeaturesMap);

        final Map<String,Map<String,Double>> cancerFeatureCounts = Maps.newHashMap();

        final Map<String,List<Double>> driversPerSampleMap = Maps.newHashMap();
        final List<Double> panCancerDriversPerSample = Lists.newArrayList();

        assignFeatures(sampleFeaturesMap, cancerFeatureCounts, driversPerSampleMap, panCancerDriversPerSample);

        writeFeaturePrevalenceFile(cancerFeatureCounts);

        writeAverageDriversFile(driversPerSampleMap, panCancerDriversPerSample);

        return true;
    }

    private boolean loadSampleData(final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        if(!mConfig.CohortFeaturesFile.isEmpty())
        {
            final Map<String,List<SampleFeatureData>> allSampleFeatures = Maps.newHashMap();

            List<String> files = parseFileSet(mConfig.CohortFeaturesFile);

            for(String file : files)
            {
                if(!loadFeaturesFromCohortFile(file, allSampleFeatures))
                    return false;
            }

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
            if(!loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.refSampleIds(false), sampleFeaturesMap))
                return false;
        }
        else
        {
            // load from per-sample files
            for(int i = 0; i < mSampleDataCache.RefSampleDataList.size(); ++i)
            {
                SampleData sample = mSampleDataCache.RefSampleDataList.get(i);

                final String linxDataDir = formSamplePath(mConfig.LinxDir, sample.Id);
                final String purpleDataDir = formSamplePath(mConfig.PurpleDir, sample.Id);
                final String virusDataDir = formSamplePath(mConfig.VirusDir, sample.Id);

                if(!loadFeaturesFromFile(sample.Id, linxDataDir, purpleDataDir, virusDataDir, sampleFeaturesMap))
                    return false;

                if(i > 0 && (i % 100) == 0)
                {
                    CUP_LOGGER.debug("processed {} sample feature files", i);
                }
            }
        }

        return true;
    }

    private void assignFeatures(
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap, final Map<String,Map<String,Double>> cancerFeatureCounts,
            final Map<String,List<Double>> driversPerSampleMap, final List<Double> panCancerDriversPerSample)
    {
        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            final String sampleId = sampleEntry.getKey();
            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

            if(!isKnownCancerType(cancerType))
                continue;

            final Set<String> sampleFeatures = Sets.newHashSet();
            double driverTotal = 0;

            for(SampleFeatureData feature : sampleEntry.getValue())
            {
                if(feature.Likelihood == 0)
                    continue;

                if(feature.Type == DRIVER || feature.Type == AMP)
                    driverTotal += feature.Likelihood;

                String nameType = featureTypeName(feature.Type, feature.Name);
                if(sampleFeatures.contains(nameType)) // keep unique contributions per sample
                    continue;

                sampleFeatures.add(nameType);

                double likelihoodTotal = sampleEntry.getValue().stream()
                        .filter(x -> x.Name.equals(feature.Name) && x.Type == feature.Type)
                        .mapToDouble(x -> x.Likelihood).sum();

                likelihoodTotal = min(likelihoodTotal, 1);

                Map<String,Double> featureCounts = cancerFeatureCounts.get(cancerType);
                if(featureCounts == null)
                {
                    featureCounts = Maps.newHashMap();
                    cancerFeatureCounts.put(cancerType, featureCounts);
                }

                Double featureCount = featureCounts.get(nameType);
                featureCounts.put(nameType, featureCount == null ? likelihoodTotal : featureCount + likelihoodTotal);
            }

            List<Double> sampleDrivers = driversPerSampleMap.get(cancerType);
            if(sampleDrivers == null)
            {
                sampleDrivers = Lists.newArrayList();
                driversPerSampleMap.put(cancerType, sampleDrivers);
            }

            sampleDrivers.add(driverTotal);
            panCancerDriversPerSample.add(driverTotal);
        }
    }

    private void writeFeaturePrevalenceFile(final Map<String,Map<String,Double>> cancerFeatureCounts)
    {
        // output: CancerType,Name,Type,SamplePerc
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

                for(Map.Entry<String,Double> featureCounts : cancerEntry.getValue().entrySet())
                {
                    String[] typeNameValues = featureCounts.getKey().split(TYPE_NAME_DELIM, 2);
                    FeatureType featureType = FeatureType.valueOf(typeNameValues[0]);
                    String featureName = typeNameValues[1];

                    // check for any feature-cancer overrides
                    final FeaturePrevData override = mFeatureOverrides.stream()
                            .filter(x -> x.CancerType.equals(cancerType))
                            .filter(x -> x.Type == featureType && x.Name.equals(featureName))
                            .findFirst().orElse(null);

                    if(override == null)
                    {
                        double sampleTotal = featureCounts.getValue();
                        double prevalence = sampleTotal / cancerSamples;

                        writer.write(String.format("%s,%s,%s,%.6f", cancerType, featureName, featureType, prevalence));
                        writer.newLine();
                    }
                    else
                    {
                        writer.write(String.format("%s,%s,%s,%.6f", cancerType, override.Name, override.Type, override.Prevalence));
                        writer.newLine();

                        // remove this entry since unobserved features will then be written at the end
                        mFeatureOverrides.remove(override);
                    }
                }
            }

            for(FeaturePrevData override : mFeatureOverrides)
            {
                writer.write(String.format("%s,%s,%s,%.6f", override.CancerType, override.Name, override.Type, override.Prevalence));
                writer.newLine();
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

        final String filename = mConfig.OutputDir + COHORT_REF_FEATURE_DATA_FILE;
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
