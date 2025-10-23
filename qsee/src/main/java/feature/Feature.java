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

    public Feature(String name, double value, @Nullable FeatureType type)
    {
        mKey = new FeatureKey(name, type);
        mValue = value;
    }

    public Feature(String name, double value)
    {
        this(name, value, null);
    }

    public FeatureKey key() { return mKey; }

    public double value() { return mValue; }
}
