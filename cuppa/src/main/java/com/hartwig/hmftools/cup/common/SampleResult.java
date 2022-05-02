package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class SampleResult
{
    public final String SampleId;
    public final CategoryType Category;
    public final ResultType Result;
    public final String DataType;
    public final String Value;

    public final Map<String,Double> CancerTypeValues;

    public SampleResult(
            final String sampleId, final CategoryType category, final ResultType resultType, final String dataType, final String value,
            final Map<String,Double> cancerTypeValues)
    {
        SampleId = sampleId;
        Category = category;
        DataType = dataType;
        Result = resultType;
        Value = value;
        CancerTypeValues = cancerTypeValues;
    }

    public static boolean checkIsValidCancerType(final SampleData sample, final String refCancerType, final Map<String,Double> cancerDataMap)
    {
        if(sample.isCandidateCancerType(refCancerType))
            return true;

        cancerDataMap.put(refCancerType, 0.0);
        return false;
    }

    public boolean typeMatches(final SampleResult other)
    {
        return SampleId.equals(other.SampleId) && Category == other.Category && DataType.equals(other.DataType)
                && Result == other.Result;
    }

    public String topRefResult()
    {
        double topValue = 0;
        String topCancerType = "";

        for(Map.Entry<String,Double> entry : CancerTypeValues.entrySet())
        {
            if(entry.getValue() > topValue)
            {
                topValue = entry.getValue();
                topCancerType = entry.getKey();
            }
        }

        return topCancerType;
    }

    public static String csvHeader()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add("SampleId");
        sj.add("Category");
        sj.add("ResultType");
        sj.add("DataType");
        sj.add("Value");
        sj.add("RefCancerType");
        sj.add("RefValue");
        return sj.toString();
    }

    public void write(final BufferedWriter writer) throws IOException
    {
        final String sampleStr = String.format("%s,%s,%s,%s,%s",
                SampleId, Category, Result, DataType, Value);

        for(Map.Entry<String,Double> cancerValues : CancerTypeValues.entrySet())
        {
            writer.write(String.format("%s,%s,%.3g",
                    sampleStr, cancerValues.getKey(), cancerValues.getValue()));
            writer.newLine();
        }
    }

    public String toString()
    {
        return String.format("sample(%s) cat(%s) resultType(%s) type(%s) value(%s)",
                SampleId, Category, Result, DataType, Value);
    }

    public static Map<String,List<SampleResult>> loadResults(final String filename)
    {
        Map<String,List<SampleResult>> allSampleResults = Maps.newHashMap();

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            final Map<String,Integer> fieldsMap  = createFieldsIndexMap(line, DATA_DELIM);

            int sampleIdIndex = fieldsMap.get("SampleId");
            int categoryIndex = fieldsMap.get("Category");
            int resultTypeIndex = fieldsMap.get("ResultType");
            int dataTypeIndex = fieldsMap.get("DataType");
            int valueIndex = fieldsMap.get("Value");
            int refCancerTypeIndex = fieldsMap.get("RefCancerType");
            int refValueIndex = fieldsMap.get("RefValue");

            String currentSampleId = "";
            List<SampleResult> sampleResults = null;
            SampleResult currentResult = null;

            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(DATA_DELIM, -1);

                String sampleId = values[sampleIdIndex];
                CategoryType category = CategoryType.valueOf(values[categoryIndex]);
                ResultType resultType = ResultType.valueOf(values[resultTypeIndex]);
                String dataType = values[dataTypeIndex];
                String value = values[valueIndex];
                String refCancerType = values[refCancerTypeIndex];
                double refValue = Double.parseDouble(values[refValueIndex]);

                if(!sampleId.equals(currentSampleId))
                {
                    sampleResults = Lists.newArrayList();
                    currentSampleId = sampleId;
                    allSampleResults.put(sampleId, sampleResults);
                    currentResult = null;
                }

                if(currentResult == null || currentResult.CancerTypeValues.containsKey(refCancerType))
                {
                    currentResult = new SampleResult(sampleId, category, resultType, dataType, value, Maps.newHashMap());
                    sampleResults.add(currentResult);
                }

                currentResult.CancerTypeValues.put(refCancerType, refValue);
            }

            CUP_LOGGER.info("loaded {} Cuppa results from samples from file({})", allSampleResults.size(), filename);
        }
        catch (IOException e)
        {
            CUP_LOGGER.warn("failed to load Cuppa results from file({}): {}", filename, e.toString());
        }

        return allSampleResults;
    }
}
