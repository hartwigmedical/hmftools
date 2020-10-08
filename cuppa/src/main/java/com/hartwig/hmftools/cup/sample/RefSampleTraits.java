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
import static com.hartwig.hmftools.cup.sample.SampleTraitsDataLoader.loadTraitsFromDatabase;

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
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class RefSampleTraits
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<SampleTraitsData>> mCancerTraitsData;

    private BufferedWriter mRefDataWriter;

    public RefSampleTraits(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerTraitsData = Maps.newHashMap();
        mRefDataWriter = null;
    }

    public void buildRefDataSets()
    {
        CUP_LOGGER.info("building sample traits reference data");

        if(mConfig.RefSampleTraitsFile.isEmpty())
        {
            if(mConfig.DbAccess == null)
                return;

            final Map<String,SampleTraitsData> sampleTraitsData = Maps.newHashMap();
            loadTraitsFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, sampleTraitsData);

            sampleTraitsData. values().forEach(x -> assignSampleTraitsData(x));
        }
        else
        {
            loadRefPurityData(mConfig.RefSampleTraitsFile);
        }

        initialiseRefDataWriter();

        for(Map.Entry<String,List<SampleTraitsData>> entry : mCancerTraitsData.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<SampleTraitsData> traitsData = entry.getValue();

            final List<Double> purityValues = traitsData.stream().map(x -> x.Purity).collect(Collectors.toList());
            writeRefDataType(cancerType, SampleTraitType.PURITY, createPercentileData(purityValues));

            final List<Double> ploidyValues = traitsData.stream().map(x -> x.Ploidy).collect(Collectors.toList());
            writeRefDataType(cancerType, SampleTraitType.PLOIDY, createPercentileData(ploidyValues));

            final List<Double> msIndelTmbValues = traitsData.stream().map(x -> x.IndelsMbPerMb).collect(Collectors.toList());
            writeRefDataType(cancerType, SampleTraitType.MS_INDELS_TMB, createPercentileData(msIndelTmbValues));

            final List<Double> chordHrdValues = traitsData.stream().map(x -> x.ChordHrd).collect(Collectors.toList());
            writeRefDataType(cancerType, SampleTraitType.CHORD_HRD, createPercentileData(chordHrdValues));
        }

        closeBufferedWriter(mRefDataWriter);
    }
    
    private void initialiseRefDataWriter()
    {
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_sample_trait_percentiles.csv";
            mRefDataWriter = createBufferedWriter(filename, false);

            mRefDataWriter.write("CancerType,TraitType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mRefDataWriter.write(String.format(",Pct_%.2f", i * 0.01));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private void writeRefDataType(final String cancerType, final SampleTraitType traitType, final double[] percentileValues)
    {
        try
        {
            mRefDataWriter.write(String.format("%s,%s", cancerType, traitType));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mRefDataWriter.write(String.format(",%.6f", percentileValues[i]));
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

    private void assignSampleTraitsData(final SampleTraitsData traitsData)
    {
        final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(traitsData.SampleId);
        if(cancerType == null)
        {
            CUP_LOGGER.error("sample({}) traits missing cancer type", traitsData.SampleId);
            return;
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

    private void loadRefPurityData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SampleTraitsData traitsData = SampleTraitsData.from(fieldsIndexMap, line);
                assignSampleTraitsData(traitsData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read ref sample traits data file({}): {}", filename, e.toString());
        }
    }


}
