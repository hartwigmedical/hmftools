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
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefSampleTraits
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<SampleTraitsData>> mCancerTraitsData;

    private BufferedWriter mPercentilesWriter;
    private BufferedWriter mRatesWriter;

    public RefSampleTraits(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerTraitsData = Maps.newHashMap();
        mPercentilesWriter = null;
        mRatesWriter = null;
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

        initialiseRefDataWriters();

        for(Map.Entry<String,List<SampleTraitsData>> entry : mCancerTraitsData.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<SampleTraitsData> traitsData = entry.getValue();

            final List<Double> purityValues = traitsData.stream().map(x -> x.Purity).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.PURITY, createPercentileData(purityValues));

            final List<Double> ploidyValues = traitsData.stream().map(x -> x.Ploidy).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.PLOIDY, createPercentileData(ploidyValues));

            final List<Double> msIndelTmbValues = traitsData.stream().map(x -> x.IndelsMbPerMb).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.MS_INDELS_TMB, createPercentileData(msIndelTmbValues));

            final List<Double> chordHrdValues = traitsData.stream().map(x -> x.ChordHrd).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.CHORD_HRD, createPercentileData(chordHrdValues));

            double cancerSamples = mSampleDataCache.RefCancerSampleData.get(cancerType).size();
            int wgdCount = (int)traitsData.stream().filter(x -> x.HasWGD).count();
            int femaleCount = (int)traitsData.stream().filter(x -> x.GenderType == Gender.FEMALE).count();

            writeRatesData(cancerType, wgdCount/cancerSamples, femaleCount/cancerSamples);
        }

        closeBufferedWriter(mPercentilesWriter);
        closeBufferedWriter(mRatesWriter);
    }
    
    private void initialiseRefDataWriters()
    {
        try
        {
            final String percFilename = mConfig.OutputDir + "cup_ref_sample_trait_percentiles.csv";
            mPercentilesWriter = createBufferedWriter(percFilename, false);

            mPercentilesWriter.write("CancerType,TraitType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mPercentilesWriter.write(String.format(",Pct_%.2f", i * 0.01));
            }

            mPercentilesWriter.newLine();

            final String ratesFilename = mConfig.OutputDir + "cup_ref_sample_trait_rates.csv";
            mRatesWriter = createBufferedWriter(ratesFilename, false);

            mRatesWriter.write("CancerType,WGDPerc,GenderFemalePerc");
            mRatesWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private void writeRatesData(final String cancerType, double wgdRate, double femaleRate)
    {
        try
        {
            mRatesWriter.write(String.format("%s,%.4f,%.4f", cancerType, wgdRate, femaleRate));
            mRatesWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref sample traits rates ref data output: {}", e.toString());
        }
    }

    private void writePercentilesData(final String cancerType, final SampleTraitType traitType, final double[] percentileValues)
    {
        try
        {
            mPercentilesWriter.write(String.format("%s,%s", cancerType, traitType));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mPercentilesWriter.write(String.format(",%.6f", percentileValues[i]));
            }

            mPercentilesWriter.newLine();
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
