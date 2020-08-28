package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.sigs.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefSvAnnotation
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<SvData>> mCancerSvData;

    public RefSvAnnotation(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerSvData = Maps.newHashMap();
    }

    public void buildRefDataSets()
    {
        if(mConfig.RefSampleSvDataFile.isEmpty())
            return;

        loadRefSvData(mConfig.RefSampleSvDataFile);

        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_sv_percentiles.csv";
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
                final List<SvData> svDataList = entry.getValue();

                for(SvDataType dataType : SvDataType.values())
                {
                    final List<Double> values = svDataList.stream().map(x -> (double)x.getCount(dataType)).collect(Collectors.toList());
                    // writeRefDataType(cancerType, dataType, createPercentileData(values));

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

    private void loadRefSvData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SvData svData = SvData.from(fieldsIndexMap, line);

                final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(svData.SampleId);
                if(cancerType == null)
                {
                    CUP_LOGGER.error("sample({}) missing cancer type", svData.SampleId);
                    continue;
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
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read ref SV data file({}): {}", filename, e.toString());
        }
    }


}
