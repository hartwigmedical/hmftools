package com.hartwig.hmftools.qsee.feature;

import static com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder.FIELD_KEY_VALUE_SEPARATOR;

import java.util.Comparator;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.prep.category.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.prep.category.discordant.DiscordantFragGroup;
import com.hartwig.hmftools.qsee.prep.category.msindel.RepeatType;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

import org.apache.commons.lang3.tuple.Pair;
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
            return typeDiff;

        int sourceToolDiff = Comparator.nullsLast(Comparator.<SourceTool>naturalOrder()).compare(this.mSourceTool, other.mSourceTool);
        if(sourceToolDiff != 0)
            return sourceToolDiff;

        if(this.mType == FeatureType.SUMMARY_TABLE)
            return SummaryTableFeature.valueOf(this.mName).compareTo(SummaryTableFeature.valueOf(other.mName));

        if(this.mType == FeatureType.DISCORDANT_FRAG_FREQ)
            return compareNameDiscordantFragFreq(this.mName, other.mName);

        if(this.mType == FeatureType.BQR_PER_ORIG_QUAL || this.mType == FeatureType.BQR_PER_SNV96_CONTEXT)
            return compareNameBqr(this.mName, other.mName);

        if(this.mType == FeatureType.MS_INDEL_ERROR_RATES || this.mType == FeatureType.MS_INDEL_ERROR_BIAS)
            return compareNameMsIndels(this.mName, other.mName);

        if(this.mName.contains(FIELD_KEY_VALUE_SEPARATOR) && other.mName.contains(FIELD_KEY_VALUE_SEPARATOR))
            return compareNameWithNumericalValues(this.mName, other.mName);

        return Comparator.nullsLast(Comparator.<String>naturalOrder()).compare(this.mName, other.mName);
    }

    private static int compareNameDiscordantFragFreq(String name, String nameOther)
    {
        DiscordantFragGroup group = DiscordantFragGroup.valueOf(name.split(FIELD_KEY_VALUE_SEPARATOR)[1]);
        DiscordantFragGroup groupOther = DiscordantFragGroup.valueOf(nameOther.split(FIELD_KEY_VALUE_SEPARATOR)[1]);

        return group.compareTo(groupOther);
    }

    private static int compareNameBqr(String name, String nameOther)
    {
        List<Pair<String, String>> nameParts = MultiFieldStringBuilder.parseMultiFieldAsPairs(name);
        List<Pair<String, String>> namePartsOther = MultiFieldStringBuilder.parseMultiFieldAsPairs(nameOther);

        BaseQualBin qualBin = BaseQualBin.valueOf(nameParts
                .get(nameParts.size()-1)
                .getValue()
                .split(" ")[0]
        );

        BaseQualBin qualBinOther = BaseQualBin.valueOf(namePartsOther
                .get(namePartsOther.size()-1)
                .getValue()
                .split(" ")[0]
        );

        return qualBin.compareTo(qualBinOther);
    }

    private static int compareNameMsIndels(String name, String nameOther)
    {
        List<Pair<String, String>> nameParts = MultiFieldStringBuilder.parseMultiFieldAsPairs(name);
        List<Pair<String, String>> namePartsOther = MultiFieldStringBuilder.parseMultiFieldAsPairs(nameOther);

        ConsensusType consensusType = ConsensusType.valueOf(nameParts.get(0).getValue());
        ConsensusType consensusTypeOther = ConsensusType.valueOf(namePartsOther.get(0).getValue());

        int consensusTypeDiff = consensusType.compareTo(consensusTypeOther);
        if(consensusTypeDiff != 0)
            return consensusTypeDiff;

        RepeatType repeatType = RepeatType.fromDisplayName(nameParts.get(1).getValue());
        RepeatType repeatTypeOther = RepeatType.fromDisplayName(namePartsOther.get(1).getValue());

        int repeatTypeDiff = repeatType.compareTo(repeatTypeOther);
        if(repeatTypeDiff != 0)
            return repeatTypeDiff;

        double numUnits = Double.parseDouble(nameParts.get(2).getValue());
        double numUnitsOther = Double.parseDouble(namePartsOther.get(2).getValue());
        return Double.compare(numUnits, numUnitsOther);
    }

    private static int compareNameWithNumericalValues(String name, String nameOther)
    {
        String[] nameParts = MultiFieldStringBuilder.parseMultiFieldAsArray(name);
        String[] namePartsOther = MultiFieldStringBuilder.parseMultiFieldAsArray(nameOther);

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
