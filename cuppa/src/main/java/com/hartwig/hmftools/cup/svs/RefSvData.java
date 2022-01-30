package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_SV_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SV_PERC;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.parseFileSet;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromCohortFile;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromDatabase;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromFile;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.ref.RefClassifier;

public class RefSvData implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<SvData>> mCancerSvData;

    public RefSvData(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerSvData = Maps.newHashMap();
    }

    public CategoryType categoryType() { return SV; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.DbAccess != null || !config.CohortSampleSvDataFile.isEmpty() || !config.SampleSvVcf.isEmpty();
    }

    public void buildRefDataSets()
    {
        if(mConfig.CohortSampleSvDataFile.isEmpty() && mConfig.SampleSvVcf.isEmpty() && mConfig.DbAccess == null)
            return;

        CUP_LOGGER.info("building SV reference data");

        final Map<String,SvData> sampleSvData = Maps.newHashMap();

        if(!mConfig.CohortSampleSvDataFile.isEmpty())
        {
            final List<String> files = parseFileSet(mConfig.CohortSampleSvDataFile);
            files.forEach(x -> loadCohortSvData(x, sampleSvData));
        }
        else if(!mConfig.SampleSvVcf.isEmpty())
        {
            // load from per-sample files
            for(int i = 0; i < mSampleDataCache.RefSampleDataList.size(); ++i)
            {
                SampleData sample = mSampleDataCache.RefSampleDataList.get(i);

                final String svVcfFile = formSamplePath(mConfig.SampleSvVcf, sample.Id);
                final String sampleDataDir = formSamplePath(mConfig.SampleFeaturesDir, sample.Id);
                final String clusterFile = LinxCluster.generateFilename(sampleDataDir, sample.Id);

                if(!loadSvDataFromFile(sample.Id, svVcfFile, clusterFile, sampleSvData))
                    break;

                if(i > 0 && (i % 100) == 0)
                {
                    CUP_LOGGER.debug("processed {} sample SV files", i);
                }
            }
        }
        else
        {
            loadSvDataFromDatabase(mConfig.DbAccess, mSampleDataCache.refSampleIds(false), sampleSvData);
            sampleSvData.values().forEach(x -> assignSampleData(x));
        }

        writeCohortData(sampleSvData);

        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_SV_PERC;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("CancerType,SvDataType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                writer.write(String.format(",Pct_%.2f", i * 0.01));
            }

            writer.newLine();

            for(Map.Entry<String,List<SvData>> entry : mCancerSvData.entrySet())
            {
                final String cancerType = entry.getKey();

                if(!isKnownCancerType(cancerType))
                    continue;

                final List<SvData> svDataList = entry.getValue();

                for(SvDataType dataType : SvDataType.values())
                {
                    final List<Double> values = svDataList.stream().map(x -> (double)x.getCount(dataType)).collect(Collectors.toList());

                    writer.write(String.format("%s,%s", cancerType, dataType));

                    final double[] percentileValues = createPercentileData(values);

                    for(int i = 0; i < percentileValues.length; ++i)
                    {
                        writer.write(String.format(",%.6f", percentileValues[i]));
                    }

                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref sample SV data output: {}", e.toString());
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

    private void assignSampleData(final SvData svData)
    {
        final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(svData.SampleId);
        if(cancerType == null)
        {
            CUP_LOGGER.error("sample({}) SV missing cancer type", svData.SampleId);
            return;
        }

        if(!isKnownCancerType(cancerType))
            return;

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

    private void writeCohortData(final Map<String,SvData> sampleSvData)
    {
        if(!mConfig.WriteCohortFiles)
            return;

        final String filename = mConfig.OutputDir + COHORT_REF_FILE_SV_DATA;

        if(Files.exists(Paths.get(filename)))
        {
            CUP_LOGGER.warn("not over-writing cohort SV reference file({})", filename);
            return;
        }

        CUP_LOGGER.info("writing cohort SV reference data");

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(SvData.header());
            writer.newLine();

            for(Map.Entry<String,SvData> entry : sampleSvData.entrySet())
            {
                final String sampleId = entry.getKey();
                final SvData svData = entry.getValue();
                writer.write(String.format("%s,%s", sampleId, svData.toCsv()));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write SV cohort data output: {}", e.toString());
        }
    }

    private void loadCohortSvData(final String filename, final Map<String,SvData> sampleSvDataMap)
    {
        final Map<String,SvData> svDataMap = Maps.newHashMap();

        loadSvDataFromCohortFile(filename, svDataMap);

        // restrict to the loaded ref samples
        for(SvData svData : svDataMap.values())
        {
            if(mSampleDataCache.hasRefSample(svData.SampleId))
            {
                sampleSvDataMap.put(svData.SampleId, svData);
                assignSampleData(svData);
            }
        }
    }

}
