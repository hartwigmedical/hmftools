package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.sigs.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefFeatures
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    public RefFeatures(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;
    }

    public void buildRefDataSets()
    {
        if(mConfig.DbAccess == null)
            return;

        final Map<String,List<SampleFeatureData>> sampleFeaturesMap = Maps.newHashMap();
        loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, sampleFeaturesMap);

        final Map<String,Map<String,Double>> cancerFeatureCounts = Maps.newHashMap();
        final Map<String,FeatureType> featureTypes = Maps.newHashMap();

        final Map<String,List<Double>> driversPerSampleMap = Maps.newHashMap();
        final List<Double> panCancerDriversPerSample = Lists.newArrayList();

        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            final String sampleId = sampleEntry.getKey();
            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

            final Set<String> features = Sets.newHashSet();
            double driverTotal = 0;

            for(SampleFeatureData feature : sampleEntry.getValue())
            {
                if(feature.Type == DRIVER)
                    driverTotal += feature.Likelihood;

                if(features.contains(feature.Name))
                    continue;

                features.add(feature.Name);
                featureTypes.put(feature.Name, feature.Type);

                double likelihoodTotal = sampleEntry.getValue().stream().filter(x -> x.Name.equals(feature.Name)).mapToDouble(x -> x.Likelihood).sum();

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

        writeFeaturePrevalenceFile(cancerFeatureCounts, featureTypes);
        writeAverageDriversFile(driversPerSampleMap, panCancerDriversPerSample);
    }

    private void writeFeaturePrevalenceFile(final Map<String,Map<String,Double>> cancerFeatureCounts, final Map<String,FeatureType> featureTypes)
    {
        // output: CancerType,Gene,Type,SamplePerc
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_feature_prev.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("CancerType,Feature,Type,SamplePerc");
            writer.newLine();

            for(Map.Entry<String,Map<String,Double>> cancerEntry : cancerFeatureCounts.entrySet())
            {
                final String cancerType = cancerEntry.getKey();

                for(Map.Entry<String,Double> featureEntry : cancerEntry.getValue().entrySet())
                {
                    final String feature = featureEntry.getKey();
                    double sampleTotal = featureEntry.getValue();
                    final FeatureType featureType = featureTypes.get(feature);

                    writer.write(String.format("%s,%s,%s,%.6f", cancerType, feature, featureType, sampleTotal));
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
            final String filename = mConfig.OutputDir + "cup_ref_avg_drivers.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("CancerType,SvDataType");
            writer.newLine();

            for(Map.Entry<String,List<Double>> entry : driversPerSampleMap.entrySet())
            {
                final String cancerType = entry.getKey();
                final List<Double> drivers = entry.getValue();

                double driverTotal = drivers.stream().mapToDouble(x -> x).sum();
                int cancerSamples = mSampleDataCache.RefCancerSampleData.get(cancerType).size();
                double avgDrivers = driverTotal / cancerSamples;

                writer.write(String.format("%s,%.2f", cancerType, avgDrivers));
                writer.newLine();
            }

            double driverTotal = panCancerDriversPerSample.stream().mapToDouble(x -> x).sum();
            double avgDrivers = driverTotal / panCancerDriversPerSample.size();

            writer.write(String.format("ALL,%.2f", avgDrivers));
            writer.newLine();

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref driver averages file: {}", e.toString());
        }
    }

    private double[] createPercentileData(final List<Double> values)
    {
        // sort the data into an array
        final List<Integer> sortedIndices = getSortedVectorIndices(convertList(values), true);
        final double[] sortedValues = new double[values.size()];

        for(int i = 0; i < values.size(); ++i)
        {
            sortedValues[i] = values.get(sortedIndices.get(i));
        }

        return buildPercentiles(sortedValues);
    }

    /*
    private void assignSampleData(final SvData svData)
    {
        final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(svData.SampleId);
        if(cancerType == null)
        {
            CUP_LOGGER.error("sample({}) missing cancer type", svData.SampleId);
            return;
        }

        List<SvData> svDataList = mCancerSvData.get(cancerType);
        if(svDataList == null)
        {
            mCancerSvData.put(cancerType, Lists.newArrayList(svData));
        }
        else
        {
            svDataList.add(svData);
        }
    }

     */
}
