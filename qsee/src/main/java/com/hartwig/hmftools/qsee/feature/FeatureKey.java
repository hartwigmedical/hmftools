package com.hartwig.hmftools.qsee.feature;

import java.util.Objects;
import java.util.StringJoiner;

public class FeatureKey
{
    private final String mName;
    private final FeatureType mType;
    private final SourceTool mSourceTool;

    // For multi-field keys
    private static final String FIELD_KEY_VALUE_SEPARATOR = "=";
    private static final String FIELD_SEPARATOR = ";";

    public FeatureKey(String name, FeatureType type, SourceTool sourceTool)
    {
        mName = name;
        mType = type;
        mSourceTool = sourceTool;
    }

    public String name() { return mName; }
    public FeatureType type() { return mType; }
    public SourceTool sourceTool() { return mSourceTool; }

    public static String formSingleFieldName(String fieldName, String fieldValue)
    {
        return fieldName + FIELD_KEY_VALUE_SEPARATOR + fieldValue;
    }

    public static String formMultiFieldName(String... keyValuePairs)
    {
        if(keyValuePairs.length % 2 != 0)
        {
            throw new IllegalArgumentException("Must provide an even number of arguments (key-value pairs)");
        }

        StringJoiner featureName = new StringJoiner(FIELD_SEPARATOR);

        for(int i = 0; i < keyValuePairs.length; i += 2)
        {
            String fieldName = keyValuePairs[i];
            String fieldValue = keyValuePairs[i + 1];
            String fieldString = fieldName + FIELD_KEY_VALUE_SEPARATOR + fieldValue;
            featureName.add(fieldString);
        }

        return featureName.toString();
    }

    @Override
    public String toString()
    {
        return (mType == null)
                ? mName
                : String.format("name(%s) type(%s)", mName, mType.name());
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        final FeatureKey other = (FeatureKey) o;
        return mName.equals(other.mName) && mType == other.mType;
    }

    @Override
    public int hashCode()
    {
        int result = Objects.hashCode(mName);
        result = 31 * result + Objects.hashCode(mType);
        return result;
    }
}
