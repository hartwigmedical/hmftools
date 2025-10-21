package feature;

import java.util.StringJoiner;

import org.apache.commons.lang3.tuple.Pair;

public class FeatureValue
{
    public String mKey;
    public double mValue;
    public FeatureType mType;

    // For multi-field keys
    private static final String KEY_VALUE_SEPARATOR = "=";
    private static final String KEY_VALUE_PAIR_SEPARATOR = ";";

    public FeatureValue(String key, double value, FeatureType type)
    {
        mKey = key;
        mValue = value;
        mType = type;
    }

    public static String keyFromPair(String fieldName, String fieldValue)
    {
        return fieldName + KEY_VALUE_SEPARATOR + fieldValue;
    }

    @SafeVarargs
    public static String keyFromPairs(Pair<String, String>... pairs)
    {
        StringJoiner joiner = new StringJoiner(KEY_VALUE_PAIR_SEPARATOR);

        for(Pair<String, String> pair : pairs)
        {
            joiner.add(keyFromPair(pair.getKey(), pair.getValue()));
        }

        return joiner.toString();
    }
}
