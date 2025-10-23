package feature;

import java.util.List;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.stream.Stream;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class FeatureKey
{
    private final String mName;
    @Nullable private final FeatureType mType;

    // For multi-field keys
    private static final String KEY_VALUE_SEPARATOR = "=";
    private static final String KEY_VALUE_PAIR_SEPARATOR = ";";

    public FeatureKey(String name, @Nullable FeatureType type)
    {
        mName = name;
        mType = type;
    }

    public FeatureKey(String name)
    {
        mName = name;
        mType = null;
    }

    public String name() { return mName; }
    public FeatureType type() { return mType; }

    public static FeatureKey of(String name) { return new FeatureKey(name); }

    public static List<FeatureKey> of(String... names) { return Stream.of(names).map(FeatureKey::new).toList(); }

    public static FeatureKey of(FeatureType type, Pair<String, String> pair)
    {
        return new FeatureKey(nameFromPair(pair.getKey(), pair.getValue()), type);
    }

    @SafeVarargs
    public static FeatureKey of(FeatureType type, Pair<String, String>... pairs)
    {
        return new FeatureKey(nameFromPairs(pairs), type);
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

    @Override
    public String toString()
    {
        return String.format("name(%s) type(%s)",
                mName,
                (mType == null) ? "UNKNOWN" : mType.toString()
        );
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
