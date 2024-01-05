package com.hartwig.hmftools.common.cuppa2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.*;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

public class CuppaPredictions
{
    public final List<CuppaPredictionEntry> PredictionEntries;
    public final List<String> CancerTypes;
    public final List<Categories.ClfName> ClfNames;
    public final List<Categories.DataType> DataTypes;
    public final Categories.ClfName MainCombinedClfName;

    public CuppaPredictions(final List<CuppaPredictionEntry> predictionEntries)
    {
        PredictionEntries = predictionEntries;

        CancerTypes = PredictionEntries.stream().map(o -> o.CancerType).distinct().collect(Collectors.toList());
        ClfNames = PredictionEntries.stream().map(o -> o.ClfName).distinct().collect(Collectors.toList());
        DataTypes = PredictionEntries.stream().map(o -> o.DataType).distinct().collect(Collectors.toList());

        MainCombinedClfName = hasRnaPredictions() ? Categories.ClfName.COMBINED : Categories.ClfName.DNA_COMBINED;
    }

    private static double parseDouble(String string)
    {
        if(string.length() == 0)
        {
            return Double.NaN;
        }

        if(string.equals("inf"))
        {
            return Double.POSITIVE_INFINITY;
        }

        if(string.equals("-inf"))
        {
            return Double.NEGATIVE_INFINITY;
        }

        return Double.parseDouble(string);
    }

    private static String parseStringWithEmpty(final String string)
    {
        if(string.length() > 0)
        {
            return string;
        }
        return "NONE";
    }

    public static CuppaPredictions fromTsv(final String filename) throws IOException
    {
        String delimiter = TSV_DELIM;

        BufferedReader fileReader = new BufferedReader(new FileReader(filename));

        String line = fileReader.readLine();
        final Map<String, Integer> fieldsMap = createFieldsIndexMap(line, delimiter);
        int sampleIdIndex = fieldsMap.get(CuppaPredictionEntry.FLD_SAMPLE_ID);
        int dataTypeIndex = fieldsMap.get(CuppaPredictionEntry.FLD_DATA_TYPE);
        int clfGroupIndex = fieldsMap.get(CuppaPredictionEntry.FLD_CLF_GROUP);
        int clfNameIndex = fieldsMap.get(CuppaPredictionEntry.FLD_CLF_NAME);
        int featNameIndex = fieldsMap.get(CuppaPredictionEntry.FLD_FEAT_NAME);
        int featValueIndex = fieldsMap.get(CuppaPredictionEntry.FLD_FEAT_VALUE);
        int cancerTypeIndex = fieldsMap.get(CuppaPredictionEntry.FLD_CANCER_TYPE);
        int dataValueIndex = fieldsMap.get(CuppaPredictionEntry.FLD_DATA_VALUE);
        int rankIndex = fieldsMap.get(CuppaPredictionEntry.FLD_RANK);
        int rankGroupIndex = fieldsMap.get(CuppaPredictionEntry.FLD_RANK_GROUP);

        List<CuppaPredictionEntry> cuppaPredictions = new ArrayList<>();
        while((line = fileReader.readLine()) != null)
        {
            String[] rowValues = line.split(delimiter, -1);

            String dataTypeStr = parseStringWithEmpty(rowValues[dataTypeIndex]).toUpperCase();
            Categories.DataType dataType = Categories.DataType.valueOf(dataTypeStr);
            if(!Categories.DataType.isSampleLevelDataType(dataType))
            {
                continue;
            }

            String sampleId = rowValues[sampleIdIndex];

            String clfGroupStr = parseStringWithEmpty(rowValues[clfGroupIndex]).toUpperCase();
            Categories.ClfGroup clfGroup = Categories.ClfGroup.valueOf(clfGroupStr);

            String clfNameStr;
            clfNameStr = parseStringWithEmpty(rowValues[clfNameIndex]).toUpperCase();
            clfNameStr = Categories.ClfName.convertAliasToName(clfNameStr);
            Categories.ClfName clfName = Categories.ClfName.valueOf(clfNameStr);

            String featName = parseStringWithEmpty(rowValues[featNameIndex]);
            double featValue = parseDouble(rowValues[featValueIndex]);
            String cancerType = rowValues[cancerTypeIndex];
            double dataValue = parseDouble(rowValues[dataValueIndex]);
            int rank = Integer.parseInt(rowValues[rankIndex]);
            int rankGroup = Integer.parseInt(rowValues[rankGroupIndex]);

            CuppaPredictionEntry cuppaPrediction = new CuppaPredictionEntry(
                    sampleId, dataType, clfGroup, clfName,
                    featName, featValue, cancerType, dataValue,
                    rank, rankGroup
            );

            cuppaPredictions.add(cuppaPrediction);
        }

        return new CuppaPredictions(cuppaPredictions);
    }

    public void printPredictions(int nRows)
    {
        int i = 0;
        for(CuppaPredictionEntry cuppaPrediction : PredictionEntries)
        {
            System.out.println( cuppaPrediction.toString());

            i++;
            if(nRows == i)
            {
                break;
            }

        }
    }

    public void printPredictions()
    {
        printPredictions(PredictionEntries.size());
    }

    public CuppaPredictionEntry get(int index)
    {
        return PredictionEntries.get(index);
    }

    public int size()
    {
        return PredictionEntries.size();
    }

    private boolean checkIsOneSample()
    {
        String targetSampleId = PredictionEntries.get(0).SampleId;

        for(CuppaPredictionEntry cuppaPrediction : PredictionEntries)
        {
            if(!cuppaPrediction.SampleId.equals(targetSampleId))
            {
                return false;
            }
        }

        return true;
    }

    private boolean hasRnaPredictions()
    {
        for(CuppaPredictionEntry cuppaPrediction : PredictionEntries)
        {
            if(!cuppaPrediction.DataType.equals(Categories.DataType.PROB))
            {
                continue;
            }

            if(cuppaPrediction.ClfName.equals(Categories.ClfName.RNA_COMBINED) & !Double.isNaN(cuppaPrediction.DataValue))
            {
                return true;
            }
        }

        return false;
    }

    public CuppaPredictions subsetByDataType(Categories.DataType dataType)
    {
        List<CuppaPredictionEntry> entries = PredictionEntries.stream()
                .filter(o -> o.DataType == dataType)
                .collect(Collectors.toList());

        return new CuppaPredictions(entries);
    }

    public CuppaPredictions subsetByClfName(Categories.ClfName clfName)
    {
        List<CuppaPredictionEntry> entries = PredictionEntries.stream()
                .filter(o -> o.ClfName == clfName)
                .collect(Collectors.toList());

        return new CuppaPredictions(entries);
    }

    public CuppaPredictions subsetByCancerType(String cancerType)
    {
        List<CuppaPredictionEntry> entries = PredictionEntries.stream()
                .filter(o -> o.CancerType.equals(cancerType))
                .collect(Collectors.toList());

        return new CuppaPredictions(entries);
    }

    public CuppaPredictions getTopPredictions(int n)
    {
        List<CuppaPredictionEntry> entries = PredictionEntries.stream()
                .filter(o -> o.Rank <= n)
                .collect(Collectors.toList());

        return new CuppaPredictions(entries);
    }

    public CuppaPredictions sortByRank()
    {
        Comparator<CuppaPredictionEntry> comparator = Comparator.comparing(cuppaPrediction -> cuppaPrediction.RankGroup);
        comparator = comparator.thenComparing(cuppaPrediction -> cuppaPrediction.Rank);

        List<CuppaPredictionEntry> sortPredictionEntries = PredictionEntries
                .stream()
                .sorted(comparator)
                .collect(Collectors.toList());

        return new CuppaPredictions(sortPredictionEntries);
    }

}


