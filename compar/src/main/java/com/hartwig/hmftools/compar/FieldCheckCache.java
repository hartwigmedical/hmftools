package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldConfigFile.FLD_ABSOLUTE_THRESHOLD;
import static com.hartwig.hmftools.compar.FieldConfigFile.FLD_COMPARED;
import static com.hartwig.hmftools.compar.FieldConfigFile.FLD_FIELD;
import static com.hartwig.hmftools.compar.FieldConfigFile.FLD_PERCENT_THRESHOLD;
import static com.hartwig.hmftools.compar.common.FileCommon.FLD_CATEGORY;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.ThresholdFieldCheck;

import org.apache.logging.log4j.Level;

public class FieldCheckCache
{
    private final Map<CategoryType,Map<String,FieldCheck>> mFieldCheckOverrides;
    private final boolean mStrictChecks;

    private final List<String> mWarnings;
    private final List<String> mErrorMessages;

    public FieldCheckCache()
    {
        this(false);
    }

    public FieldCheckCache(boolean strictChecks)
    {
        mStrictChecks = strictChecks;
        mFieldCheckOverrides = Maps.newHashMap();
        mWarnings = Lists.newArrayList();
        mErrorMessages = Lists.newArrayList();
    }

    public Map<String,FieldCheck> getCategoryFieldCheckOverrides(final CategoryType category)
    {
        return mFieldCheckOverrides.containsKey(category) ? mFieldCheckOverrides.get(category) : Collections.emptyMap();
    }

    public List<String> warnings() { return mWarnings; }
    public List<String> errorMessages() { return mErrorMessages; }

    private static final String NOT_ASSESSED_NA = "N/A";

    public void loadOverrides(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
        int fieldIndex = fieldsIndexMap.get(FLD_FIELD);
        int comparedIndex = fieldsIndexMap.get(FLD_COMPARED);
        int absoluteThresholdIndex = fieldsIndexMap.get(FLD_ABSOLUTE_THRESHOLD);
        int percentThresholdIndex = fieldsIndexMap.get(FLD_PERCENT_THRESHOLD);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            try
            {
                CategoryType categoryType = CategoryType.valueOf(values[categoryIndex]);
                String fieldName = values[fieldIndex];
                boolean isCompared = Boolean.parseBoolean(values[comparedIndex]);

                Double absThreshold = parseThreshold(values[absoluteThresholdIndex]);
                Double percThreshold = parseThreshold(values[percentThresholdIndex]);

                FieldCheck fieldCheck;

                if(absThreshold != null || percThreshold != null)
                {
                    fieldCheck = new ThresholdFieldCheck(isCompared, absThreshold, percThreshold);
                }
                else
                {
                    fieldCheck = new FieldCheck(isCompared);
                }

                Map<String,FieldCheck> categoryOverrides = mFieldCheckOverrides.get(categoryType);

                if(categoryOverrides == null)
                {
                    categoryOverrides = Maps.newHashMap();
                    mFieldCheckOverrides.put(categoryType, categoryOverrides);
                }

                categoryOverrides.put(fieldName, fieldCheck);
            }
            catch(Exception e)
            {
                Level logLevel = mStrictChecks ? Level.ERROR : Level.WARN;

                String message = format("invalid field override entry: %s", line);
                CMP_LOGGER.log(logLevel, message);

                recordProblem(message, mStrictChecks);
            }
        }

        CMP_LOGGER.info("loaded {} field check overrides from {}", lines.size(), filename);
    }

    private static Double parseThreshold(final String value)
    {
        if(value.isEmpty() || value.equals(NOT_ASSESSED_NA))
            return null;

        return Double.parseDouble(value);
    }

    public boolean hasErrors() { return !mErrorMessages.isEmpty(); }

    public void logProblems()
    {
        mWarnings.forEach(CMP_LOGGER::warn);
        mErrorMessages.forEach(CMP_LOGGER::error);
    }

    public static ThresholdFieldCheck getOrMakeFieldCheck(
            final Map<String,FieldCheck> fieldCheckMap, final String field, final Double absThreshold, final Double percThreshold)
    {
        FieldCheck override = fieldCheckMap.get(field);
        return override != null ? (ThresholdFieldCheck)override : new ThresholdFieldCheck(true, absThreshold, percThreshold);
    }

    public static FieldCheck getOrMakeFieldCheck(final Map<String,FieldCheck> fieldCheckMap, final String field)
    {
        FieldCheck override = fieldCheckMap.get(field);
        return override != null ? override : new FieldCheck(true);
    }

    public void validateAllComparisonFieldsAreSet(final ComparConfig config)
    {
        if(!mStrictChecks)
            return;

        List<ItemComparer> comparers = ComparerUtils.buildComparers(config, null);

        for(ItemComparer comparer : comparers)
        {
            List<FieldInfo> fields = comparer.fieldsList();

            Map<String,FieldCheck> categoryOverrdes = mFieldCheckOverrides.get(comparer.category());

            for(FieldInfo fieldInfo : fields)
            {
                if(!fieldInfo.fieldCheck().IsCompared) // ignore display only fields
                    continue;

                if(categoryOverrdes == null || !categoryOverrdes.containsKey(fieldInfo.name()))
                {
                    String message = format("category(%s) field(%s) missing from overriedes file", comparer.category(), fieldInfo.name());
                    recordProblem(message, true);
                }
            }
        }
    }

    private void recordProblem(final String message, final boolean isError)
    {
        if(isError)
        {
            mErrorMessages.add(message);
        }
        else
        {
            mWarnings.add(message);
        }
    }

    public void writeFieldOverridesFile(final ComparConfig config)
    {
        CMP_LOGGER.info("write field config file");
        List<ItemComparer> comparers = ComparerUtils.buildComparers(config, null);

        String filename = FieldConfigFile.generateFileName(config.formOutputFilePrefix());
        FieldConfigFile.write(filename, comparers);
    }
}
