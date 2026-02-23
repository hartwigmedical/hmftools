package com.hartwig.hmftools.qsee.feature;

import java.util.Objects;

public class FeatureKey
{
    private final String mName;
    private final FeatureType mType;
    private final SourceTool mSourceTool;

    public FeatureKey(String name, FeatureType type, SourceTool sourceTool)
    {
        mName = name;
        mType = type;
        mSourceTool = sourceTool;
    }

    public String name() { return mName; }
    public FeatureType type() { return mType; }
    public SourceTool sourceTool() { return mSourceTool; }

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
