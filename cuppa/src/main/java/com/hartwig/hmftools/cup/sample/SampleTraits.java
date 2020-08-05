package com.hartwig.hmftools.cup.sample;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.sigs.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DB_FILE_DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class SampleTraits
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;
    private final Map<String,List<SampleTraitsData>> mCancerTraitsData;

    private BufferedWriter mRefDataWriter;

    public SampleTraits(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerTraitsData = Maps.newHashMap();
        mRefDataWriter = null;

        initialiseRefDataWriter();

        loadRefPurityData(mConfig.SampleTraitsFile);
    }

    public void buildRefDataSets()
    {
        for(Map.Entry<String,List<SampleTraitsData>> entry : mCancerTraitsData.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<SampleTraitsData> traitsData = entry.getValue();

            final List<Double> purityValues = traitsData.stream().map(x -> x.Purity).collect(Collectors.toList());
            double[] percentileValues = createPercentileData(purityValues);

            writeRefDataType(cancerType, "Purity", percentileValues);
        }

        closeBufferedWriter(mRefDataWriter);
    }
    
    private void initialiseRefDataWriter()
    {
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_sample_traits.csv";
            mRefDataWriter = createBufferedWriter(filename, false);

            mRefDataWriter.write("CancerType,DataType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mRefDataWriter.write(String.format(",%.2f", i * 0.01));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private void writeRefDataType(final String cancerType, final String dataType, final double[] percentileValues)
    {
        try
        {
            mRefDataWriter.write(String.format("%s,%s", cancerType, dataType));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mRefDataWriter.write(String.format(",Pct_%.6f", percentileValues[i]));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
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

    private void loadRefPurityData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);
            int sampleIndex = fieldsIndexMap.get("SampleId");
            int genderIndex = fieldsIndexMap.get("Gender");
            int wgdIndex = fieldsIndexMap.get("WholeGenomeDuplication");
            int purityIndex = fieldsIndexMap.get("Purity");
            int ploidyIndex = fieldsIndexMap.get("Ploidy");
            int snvMbIndex = fieldsIndexMap.get("TmbPerMb");
            int msiScoreIndex = fieldsIndexMap.get("MsIndelsPerMb");

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                SampleTraitsData traitsData = new SampleTraitsData(
                        items[sampleIndex], Gender.valueOf(items[genderIndex]), items[wgdIndex].equals("1"),
                        Double.parseDouble(items[purityIndex]), Double.parseDouble(items[ploidyIndex]),
                        Double.parseDouble(items[snvMbIndex]), Double.parseDouble(items[msiScoreIndex]));

                final String cancerType = mSampleDataCache.SampleCancerTypeMap.get(traitsData.SampleId);
                if(cancerType == null)
                {
                    CUP_LOGGER.error("sample({}) missing cancer type", traitsData.SampleId);
                    continue;
                }

                List<SampleTraitsData> traitsList = mCancerTraitsData.get(cancerType);
                if(traitsList == null)
                {
                    mCancerTraitsData.put(cancerType, Lists.newArrayList(traitsData));
                }
                else
                {
                    traitsList.add(traitsData);
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read driver prevalence data file({}): {}", filename, e.toString());
        }
    }


}
