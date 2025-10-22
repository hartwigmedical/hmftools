package feature;

import java.util.StringJoiner;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class Feature
{
    public final String mKey;
    public final double mValue;

    @Nullable
    public final FeatureType mType;

    // For multi-field keys
    private static final String KEY_VALUE_SEPARATOR = "=";
    private static final String KEY_VALUE_PAIR_SEPARATOR = ";";

    public Feature(String key, double value, @Nullable FeatureType type)
    {
        mKey = key;
        mValue = value;
        mType = type;
    }

    public Feature(String key, double value)
    {
        mKey = key;
        mValue = value;
        mType = null;
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
