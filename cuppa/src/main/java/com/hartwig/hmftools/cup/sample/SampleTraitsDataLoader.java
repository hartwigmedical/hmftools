package com.hartwig.hmftools.cup.sample;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.sample.SampleTraitType.GENDER;
import static com.hartwig.hmftools.cup.sample.SampleTraitType.WGD;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SampleTraitsDataLoader
{
    public static boolean loadTraitsFromCohortFile(final String filename, final Map<String,SampleTraitsData> sampleTraitsData)
    {
        if(filename == null)
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SampleTraitsData traitsData = SampleTraitsData.from(fieldsIndexMap, line);
                sampleTraitsData.put(traitsData.SampleId, traitsData);
            }
        } catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample traits data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadTraitsFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,SampleTraitsData> sampleTraitsData)
    {
        if(dbAccess == null)
            return false;

        for(final String sampleId : sampleIds)
        {
            final PurityContext purityContext = dbAccess.readPurityContext(sampleId);
            if(purityContext == null)
            {
                CUP_LOGGER.warn("sample({}) missing purity data", sampleId);
                continue;
            }

            // ignored for now
            //final ChordAnalysis chordAnalysis = dbAccess.readChordAnalysis(sampleId);
            double chordHrd = 0; // chordAnalysis != null ? chordAnalysis.hrdValue() : 0;

            SampleTraitsData traitsData = SampleTraitsData.from(sampleId, purityContext, chordHrd);
            sampleTraitsData.put(traitsData.SampleId, traitsData);
        }

        return true;
    }

    public static boolean loadRefPercentileData(final String filename, final Map<SampleTraitType,Map<String,double[]>> refTraitPercentiles)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String cancerType = items[0];
                final SampleTraitType traitType = SampleTraitType.valueOf(items[1]);

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                Map<String,double[]> traitData = refTraitPercentiles.get(traitType);

                if(traitData == null)
                {
                    traitData = Maps.newHashMap();
                    refTraitPercentiles.put(traitType, traitData);
                }

                traitData.put(cancerType, percentileData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample traits perc data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadRefRateData(final String filename, final Map<SampleTraitType,Map<String,Double>> refTraitRates)
    {
        // CancerType,IsFemale,WGD,SampleCount,GenderFemalePerc,WGDPerc
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int wgdRateIndex = fieldsIndexMap.get("WGDPerc");
            int genderIndex = fieldsIndexMap.get("GenderFemalePerc");

            Map<String,Double> wgdRates = Maps.newHashMap();
            Map<String,Double> genderRates = Maps.newHashMap();

            refTraitRates.put(WGD, wgdRates);
            refTraitRates.put(GENDER, genderRates);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String cancerType = items[0];

                double wgdRate = Double.parseDouble(items[wgdRateIndex]);
                double genderFemaleRate = Double.parseDouble(items[genderIndex]);

                wgdRates.put(cancerType, wgdRate);
                genderRates.put(cancerType, genderFemaleRate);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample traits rate data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

}
