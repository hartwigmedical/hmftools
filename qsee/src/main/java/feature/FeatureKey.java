package feature;

import java.util.List;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class FeatureKey
{
    @Nullable private final FeatureType mType;
    private final String mName;

    // For multi-field keys
    private static final String KEY_VALUE_SEPARATOR = "=";
    private static final String KEY_VALUE_PAIR_SEPARATOR = ";";

    public FeatureKey(@Nullable FeatureType type, String name)
    {
        mType = type;
        mName = name;
    }

    public FeatureKey(String name)
    {
        mType = null;
        mName = name;
    }

    public String name() { return mName; }
    public FeatureType type() { return mType; }

    public static FeatureKey of(FeatureType type, String name) { return new FeatureKey(type, name); }

    public static FeatureKey of(String name) { return new FeatureKey(name); }

    @VisibleForTesting
    public static List<FeatureKey> ofNames(String... names) { return Stream.of(names).map(FeatureKey::new).toList(); }

    public static FeatureKey ofPair(FeatureType type, Pair<String, String> pair)
    {
        return new FeatureKey(type, nameFromPair(pair.getKey(), pair.getValue()));
    }

    @SafeVarargs
    public static FeatureKey ofPairs(FeatureType type, Pair<String, String>... pairs)
    {
        return new FeatureKey(type, nameFromPairs(pairs));
    }

    @SafeVarargs
    public static FeatureKey ofPairs(Pair<String, String>... pairs)
    {
        return new FeatureKey(null, nameFromPairs(pairs));
    }

    private static String nameFromPair(String fieldName, String fieldValue)
    {
        return fieldName + KEY_VALUE_SEPARATOR + fieldValue;
    }

    @SafeVarargs
    private static String nameFromPairs(Pair<String, String>... pairs)
    {
        StringJoiner joiner = new StringJoiner(KEY_VALUE_PAIR_SEPARATOR);

        for(Pair<String, String> pair : pairs)
        {
            String fieldString = nameFromPair(pair.getKey(), pair.getValue());
            joiner.add(fieldString);
        }

        return joiner.toString();
    }

    public FeatureKey withType(FeatureType type) { return new FeatureKey(type, mName); }

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
