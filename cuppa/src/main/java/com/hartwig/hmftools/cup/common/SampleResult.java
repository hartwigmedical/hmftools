package com.hartwig.hmftools.cup.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.COMBINED;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.ClassifierType.ALT_SJ_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENDER;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE;

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
            final String sampleId, final CategoryType category, final ResultType resultType,
            final String dataType, final String value, final Map<String,Double> cancerTypeValues)
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

    public static String detailedHeader()
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

    public void writeDetailed(final BufferedWriter writer) throws IOException
    {
        final String sampleStr = format("%s,%s,%s,%s,%s",
                SampleId, Category, Result, DataType, Value);

        for(Map.Entry<String,Double> cancerValues : CancerTypeValues.entrySet())
        {
            writer.write(format("%s,%s,%.3g",
                    sampleStr, cancerValues.getKey(), cancerValues.getValue()));
            writer.newLine();
        }
    }

    public String toString()
    {
        return format("sample(%s) cat(%s) resultType(%s) type(%s) value(%s)",
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

                String dataType = values[dataTypeIndex];

                if(dataType.equals(GENDER.toString()))
                    continue;

                String categoryStr = values[categoryIndex];
                ResultType resultType;

                CategoryType category;

                if(categoryStr.equals("CLASSIFIER"))
                {
                    resultType = ResultType.CLASSIFIER;

                    if(dataType.equals("SNV_96_PAIRWISE_SIMILARITY"))
                        dataType = SNV_96_PAIRWISE.toString();
                    else if(dataType.equals("GENOMIC_POSITION_SIMILARITY"))
                        dataType = GENOMIC_POSITION_COHORT.toString();

                    // translate old types:
                    if(dataType.equals(SNV_96_PAIRWISE.toString()))
                        category = SNV;
                    else if(dataType.equals(GENOMIC_POSITION_COHORT.toString()))
                        category = SNV;
                    else if(dataType.equals(ClassifierType.FEATURE.toString()))
                        category = FEATURE;
                    else if(dataType.equals(ALT_SJ_COHORT.toString()))
                        category = ALT_SJ;
                    else if(dataType.equals(EXPRESSION_PAIRWISE.toString()))
                        category = GENE_EXP;
                    else
                        continue;

                }
                else
                {
                    category = CategoryType.valueOf(categoryStr);
                    resultType = ResultType.valueOf(values[resultTypeIndex]);

                    if(category == COMBINED)
                        resultType = ResultType.CLASSIFIER;
                }

                if(resultType != ResultType.CLASSIFIER && category != COMBINED)
                    continue;

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
