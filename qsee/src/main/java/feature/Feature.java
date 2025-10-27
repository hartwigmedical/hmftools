package feature;

import org.jetbrains.annotations.Nullable;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
    }

    public Feature(@Nullable FeatureType type, String name, double value)
    {
        mKey = new FeatureKey(type, name);
        mValue = value;
    }

    public Feature(String name, double value)
    {
        this(null, name, value);
    }

    public FeatureKey key() { return mKey; }

    public double value() { return mValue; }
}
