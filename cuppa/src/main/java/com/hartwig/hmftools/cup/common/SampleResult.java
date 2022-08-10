package com.hartwig.hmftools.cup.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_CATEGORY;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_DATA_TYPE;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_REF_CANCER_TYPE;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_REF_VALUE;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_RESULT_TYPE;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.FLD_VALUE;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.isOldCategoryClassifierType;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.mapOldCategoryType;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.mapOldDataType;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.GENDER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.ResultType;

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

    public CuppaDataFile toCuppaData()
    {
        return new CuppaDataFile(Category, Result, DataType, Value, CancerTypeValues);
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
            int categoryIndex = fieldsMap.get(FLD_CATEGORY);
            int resultTypeIndex = fieldsMap.get(FLD_RESULT_TYPE);
            int dataTypeIndex = fieldsMap.get(FLD_DATA_TYPE);
            int valueIndex = fieldsMap.get(FLD_VALUE);
            int refCancerTypeIndex = fieldsMap.get(FLD_REF_CANCER_TYPE);
            int refValueIndex = fieldsMap.get(FLD_REF_VALUE);

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

                if(isOldCategoryClassifierType(categoryStr)) // support for pre-1.7
                {
                    dataType = mapOldDataType(dataType);
                    resultType = ResultType.CLASSIFIER;
                    category = mapOldCategoryType(dataType);

                    if(category == null)
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

            CUP_LOGGER.info("loaded {} Cuppa results for samples from file({})", allSampleResults.size(), filename);
        }
        catch (IOException e)
        {
            CUP_LOGGER.warn("failed to load Cuppa results from file({}): {}", filename, e.toString());
        }

        return allSampleResults;
    }
}
