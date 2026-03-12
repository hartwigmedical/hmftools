package com.hartwig.hmftools.qsee.feature;

import static com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder.FIELD_KEY_VALUE_SEPARATOR;

import java.util.Comparator;
import java.util.Objects;

import com.hartwig.hmftools.qsee.prep.category.discordant.DiscordantFragGroup;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

import org.jetbrains.annotations.NotNull;

public class FeatureKey implements Comparable<FeatureKey>
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

    @Override
    public int compareTo(@NotNull final FeatureKey other)
    {
        int typeDiff = Comparator.nullsLast(Comparator.<FeatureType>naturalOrder()).compare(this.mType, other.mType);
        if(typeDiff != 0)
        {
            return typeDiff;
        }

        int sourceToolDiff = Comparator.nullsLast(Comparator.<SourceTool>naturalOrder()).compare(this.mSourceTool, other.mSourceTool);
        if(sourceToolDiff != 0)
        {
            return sourceToolDiff;
        }

        // Handle feature names that are composed of enums
        if(this.mType == FeatureType.SUMMARY_TABLE)
        {
            return SummaryTableFeature.valueOf(this.mName).compareTo(SummaryTableFeature.valueOf(other.mName));
        }

        if(this.mType == FeatureType.DISCORDANT_FRAG_FREQ)
        {
            DiscordantFragGroup discordantFragGroup = DiscordantFragGroup.valueOf(this.mName.split(FIELD_KEY_VALUE_SEPARATOR)[1]);
            DiscordantFragGroup discordantFragGroupOther = DiscordantFragGroup.valueOf(other.mName.split(FIELD_KEY_VALUE_SEPARATOR)[1]);
            return discordantFragGroup.compareTo(discordantFragGroupOther);
        }

        // Default to numeric feature value aware string comparison
        boolean comparingMultiFields = this.mName.contains(FIELD_KEY_VALUE_SEPARATOR) && other.mName.contains(FIELD_KEY_VALUE_SEPARATOR);
        if(!comparingMultiFields)
        {
            return Comparator.nullsLast(Comparator.<String>naturalOrder()).compare(this.mName, other.mName);
        }

        String[] nameParts = this.mName.split(FIELD_KEY_VALUE_SEPARATOR);
        String[] namePartsOther = other.mName.split(FIELD_KEY_VALUE_SEPARATOR);

        int minNumParts = Math.min(nameParts.length, namePartsOther.length);

        for(int i = 0; i < minNumParts; i++)
        {
            String part = nameParts[i];
            String partOther = namePartsOther[i];

            try
            {
                double numericPart = Double.parseDouble(part);
                double numericPartOther = Double.parseDouble(partOther);

                int numericDiff = Double.compare(numericPart, numericPartOther);
                if(numericDiff != 0)
                {
                    return numericDiff;
                }
            }
            catch (NumberFormatException e)
            {
                int stringDiff = part.compareTo(partOther);
                if(stringDiff != 0)
                {
                    return stringDiff;
                }
            }
        }

        return Integer.compare(nameParts.length, namePartsOther.length);
    }
}
