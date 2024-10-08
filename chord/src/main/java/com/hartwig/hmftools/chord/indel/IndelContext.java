package com.hartwig.hmftools.chord.indel;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class IndelContext
{
    public MutationType mMutationType;
    public ContextType mContextType;
    public int mHomBasesCount;
    public int mIndelLength;

    public static final int INDEL_LENGTH_CAP = 5;
    public static final int HOM_BASES_COUNT_CAP = 5;

    public IndelContext(MutationType mutationType, ContextType contextType, int homBasesCount, int indelLength)
    {
        mMutationType = mutationType;
        mContextType = contextType;

        mHomBasesCount = homBasesCount;
        mIndelLength = indelLength;
    }

    public static IndelContext from(IndelDetails indelDetails)
    {
        return new IndelContext(
                indelDetails.mMutationType,
                indelDetails.mContextType,
                indelDetails.mMaxHomBasesCount,
                indelDetails.mIndelLength
        );
    }

    public static Map<String, Integer> initializeCounts()
    {
        List<MutationType> mutationTypes = List.of(MutationType.values());
        List<ContextType> contextTypes = List.of(ContextType.values());

        Map<String, Integer> contextCounts = new LinkedHashMap<>();

        for(ContextType contextType : contextTypes)
        {
            for(MutationType mutationType : mutationTypes)
            {
                int lengthTypeCap = (contextType == ContextType.MICROHOMOLOGY) ? HOM_BASES_COUNT_CAP : INDEL_LENGTH_CAP;

                for(int i = 1; i <= lengthTypeCap; i++)
                {
                    String contextName = mutationType.Abbreviation + contextType.Suffix + i;
                    contextCounts.put(contextName, 0);
                }
            }
        }

        return contextCounts;
    }

    public String getContextName()
    {
        int lenthTypeValueCapped = (mContextType == ContextType.MICROHOMOLOGY) ?
            Math.min(mHomBasesCount, HOM_BASES_COUNT_CAP) :
            Math.min(mIndelLength, INDEL_LENGTH_CAP);

        return mMutationType.Abbreviation + mContextType.Suffix + lenthTypeValueCapped;
    }
}
