package com.hartwig.hmftools.common.cuppa2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.*;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

public class CuppaVisData
{
    public final List<CuppaVisDataEntry> VisDataEntries;

    public CuppaVisData(final List<CuppaVisDataEntry> visDataEntries)
    {
        VisDataEntries = visDataEntries;
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

    public static CuppaVisData fromTsv(final String filename) throws IOException
    {
        String delimiter = TSV_DELIM;

        BufferedReader fileReader = new BufferedReader(new FileReader(filename));

        String line = fileReader.readLine();
        final Map<String, Integer> fieldsMap = createFieldsIndexMap(line, delimiter);
        int sampleIdIndex = fieldsMap.get(CuppaVisDataEntry.FLD_SAMPLE_ID);
        int dataTypeIndex = fieldsMap.get(CuppaVisDataEntry.FLD_DATA_TYPE);
        int clfGroupIndex = fieldsMap.get(CuppaVisDataEntry.FLD_CLF_GROUP);
        int clfNameIndex = fieldsMap.get(CuppaVisDataEntry.FLD_CLF_NAME);
        int featNameIndex = fieldsMap.get(CuppaVisDataEntry.FLD_FEAT_NAME);
        int featValueIndex = fieldsMap.get(CuppaVisDataEntry.FLD_FEAT_VALUE);
        int cancerTypeIndex = fieldsMap.get(CuppaVisDataEntry.FLD_CANCER_TYPE);
        int dataValueIndex = fieldsMap.get(CuppaVisDataEntry.FLD_DATA_VALUE);
        int rankIndex = fieldsMap.get(CuppaVisDataEntry.FLD_RANK);
        int rankGroupIndex = fieldsMap.get(CuppaVisDataEntry.FLD_RANK_GROUP);

        List<CuppaVisDataEntry> visDataEntries = new ArrayList<>();
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

            CuppaVisDataEntry visDataEntry = new CuppaVisDataEntry(
                    sampleId, dataType, clfGroup, clfName,
                    featName, featValue, cancerType, dataValue,
                    rank, rankGroup
            );

            visDataEntries.add(visDataEntry);
        }

        return new CuppaVisData(visDataEntries);
    }

    public void printEntries(int nRows)
    {
        int i = 0;
        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            System.out.println( visDataEntry.toString());

            i++;
            if(nRows == i)
            {
                break;
            }

        }
    }

    public void printEntries()
    {
        printEntries(10);
    }

    public CuppaVisDataEntry get(int index)
    {
        return VisDataEntries.get(index);
    }

    public int size()
    {
        return VisDataEntries.size();
    }

    private boolean checkIsOneSample()
    {
        String targetSampleId = VisDataEntries.get(0).SampleId;

        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            if(!visDataEntry.SampleId.equals(targetSampleId))
            {
                return false;
            }
        }

        return true;
    }

    public boolean checkHasRnaPredictions()
    {
        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            if(!visDataEntry.DataType.equals(Categories.DataType.PROB))
            {
                continue;
            }

            if(visDataEntry.ClfName.equals(Categories.ClfName.RNA_COMBINED) & !Double.isNaN(visDataEntry.DataValue))
            {
                return true;
            }
        }

        return false;
    }

    public Categories.ClfName getMainCombinedClfName()
    {
        if(checkHasRnaPredictions())
        {
            return Categories.ClfName.COMBINED;
        }
        return Categories.ClfName.DNA_COMBINED;
    }

    public CuppaVisData subsetByDataType(Categories.DataType dataType)
    {
        List<CuppaVisDataEntry> newVisDataEntries = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            if(!visDataEntry.DataType.equals(dataType))
            {
                continue;
            }
            newVisDataEntries.add(visDataEntry);
        }
        return new CuppaVisData(newVisDataEntries);
    }

    public CuppaVisData subsetByClfName(Categories.ClfName clfName)
    {
        List<CuppaVisDataEntry> newVisDataEntries = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            if(!visDataEntry.ClfName.equals(clfName))
            {
                continue;
            }
            newVisDataEntries.add(visDataEntry);
        }
        return new CuppaVisData(newVisDataEntries);
    }

    public CuppaVisData getTopPredictions(int n)
    {
        List<CuppaVisDataEntry> newVisDataEntries = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : VisDataEntries)
        {
            if(visDataEntry.Rank <= n)
            {
                newVisDataEntries.add(visDataEntry);
            }
        }
        return new CuppaVisData(newVisDataEntries);
    }

    public CuppaVisData sortByRank()
    {
        Comparator<CuppaVisDataEntry> comparator = Comparator.comparing(visDataEntry -> visDataEntry.RankGroup);
        comparator = comparator.thenComparing(visDataEntry -> visDataEntry.Rank);

        List<CuppaVisDataEntry> sortedVisDataEntries = VisDataEntries
                .stream()
                .sorted(comparator)
                .collect(Collectors.toList());

        return new CuppaVisData(sortedVisDataEntries);
    }

}


