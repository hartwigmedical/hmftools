package com.hartwig.hmftools.common.cuppa2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

public class CuppaDataFile
{
    public final String Filename;
    public final List<CuppaPrediction> CuppaPredictions;
    public final boolean HasRnaData;
    public final Categories.ClfName MainCombinedClf;
    public final List<Entry<String, Double>> SortedCancerTypeProbs;

    CuppaDataFile(final String filename) throws IOException
    {
        Filename = filename;
        CuppaPredictions = readTable(filename);
        HasRnaData = checkHasRnaData();
        MainCombinedClf = getMainCombinedClf();
        SortedCancerTypeProbs = sortMapByValue(getCancerTypeProbs(MainCombinedClf));
    }

    private static double parseDouble(String string)
    {
        if(string.length() == 0)
        {
            string = "NaN";
        }
        else if(string.equals("inf"))
        {
            string = "Infinity";
        }
        else if(string.equals("-inf"))
        {
            string = "-Infinity";
        }

        return Double.parseDouble(string);
    }

    private static String parseString(final String string)
    {
        if(string.length() > 0)
        {
            return string;
        }
        return "none";
    }

    private static List<CuppaPrediction> readTable(final String filename) throws IOException
    {
        String delimiter = TSV_DELIM;

        BufferedReader fileReader = new BufferedReader(new FileReader(filename));

        String line = fileReader.readLine();
        final Map<String, Integer> fieldsMap = createFieldsIndexMap(line, delimiter);
        int sampleIdIndex = fieldsMap.get(CuppaPrediction.FLD_SAMPLE_ID);
        int dataTypeIndex = fieldsMap.get(CuppaPrediction.FLD_DATA_TYPE);
        int clfGroupIndex = fieldsMap.get(CuppaPrediction.FLD_CLF_GROUP);
        int clfNameIndex = fieldsMap.get(CuppaPrediction.FLD_CLF_NAME);
        int featNameIndex = fieldsMap.get(CuppaPrediction.FLD_FEAT_NAME);
        int featValueIndex = fieldsMap.get(CuppaPrediction.FLD_FEAT_VALUE);
        int cancerTypeIndex = fieldsMap.get(CuppaPrediction.FLD_CANCER_TYPE);
        int dataValueIndex = fieldsMap.get(CuppaPrediction.FLD_DATA_VALUE);

        List<CuppaPrediction> cuppaPredictions = new ArrayList<>();
        while((line = fileReader.readLine()) != null)
        {
            String[] rowValues = line.split(delimiter, -1);

            String sampleId = parseString(rowValues[sampleIdIndex]);

            String dataTypeStr = parseString(rowValues[dataTypeIndex]).toUpperCase();
            Categories.DataType dataType = Categories.DataType.valueOf(dataTypeStr);

            String clfGroupStr = parseString(rowValues[clfGroupIndex]).toUpperCase();
            Categories.ClfGroup clfGroup = Categories.ClfGroup.valueOf(clfGroupStr);

            String clfNameStr = parseString(rowValues[clfNameIndex]).toUpperCase();
            Categories.ClfName clfName = Categories.ClfName.valueOf(clfNameStr);

            String featName = parseString(rowValues[featNameIndex]);
            double featValue = parseDouble(rowValues[featValueIndex]);
            String cancerType = parseString(rowValues[cancerTypeIndex]);
            double dataValue = parseDouble(rowValues[dataValueIndex]);

            CuppaPrediction cuppaPrediction = new CuppaPrediction(
                    sampleId, dataType, clfGroup, clfName,
                    featName, featValue, cancerType, dataValue
            );

            cuppaPredictions.add(cuppaPrediction);
        }

        return cuppaPredictions;
    }

    private boolean checkHasRnaData()
    {
        for(CuppaPrediction cuppaPrediction : CuppaPredictions)
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

    private Categories.ClfName getMainCombinedClf()
    {
        if(HasRnaData)
        {
            return Categories.ClfName.COMBINED;
        }
        return Categories.ClfName.DNA_COMBINED;
    }

    public LinkedHashMap<String, Double> getCancerTypeProbs(Categories.ClfName targetClfName)
    {
        LinkedHashMap<String, Double> cancerTypeProbs = new LinkedHashMap<>();

        for(CuppaPrediction cuppaPrediction : CuppaPredictions)
        {
            boolean isProbDataType = cuppaPrediction.DataType.equals(Categories.DataType.PROB);
            boolean isTargetClfName = cuppaPrediction.ClfName.equals(targetClfName);

            if(!isProbDataType | !isTargetClfName)
            {
                continue;
            }

            cancerTypeProbs.put(cuppaPrediction.CancerType, cuppaPrediction.DataValue);
        }

        return cancerTypeProbs;
    }

    public static List<Entry<String, Double>> sortMapByValue(LinkedHashMap<String, Double> map)
    {
        List<Map.Entry<String, Double>> list = new ArrayList<>(map.entrySet());
        list.sort(Map.Entry.<String, Double>comparingByValue().reversed());

        List<Entry<String, Double>> sortedList = new ArrayList<>(list);
        return sortedList;
    }
}


