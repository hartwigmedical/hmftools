package com.hartwig.hmftools.geneutils.fusion;

import static java.lang.String.format;

import com.hartwig.hmftools.common.fusion.KnownFusionType;

public class FusionRefData
{
    public final KnownFusionType Type;
    public final String FiveGene;
    public final String ThreeGene;
    public final String CancerTypes;
    public final String PubMedId;
    public final String KnownExonTranscript37;
    public final String KnownExonUpRange37;
    public final String KnownExonDownRange37;
    public final String HighImpactPromiscuous;
    public final String OverridesRef37;
    public final String KnownExonTranscriptRef38;
    public final String KnownExonUpRangeRef38;
    public final String KnownExonDownRangeRef38;
    public final String OverridesRef38;

    public FusionRefData(
            final KnownFusionType type, final String fiveGene, final String threeGene, final String cancerTypes, final String pubMedId,
            final String knownExonTranscript37, final String knownExonUpRange37, final String knownExonDownRange37,
            final String highImpactPromiscuous, final String overridesRef37, final String knownExonTranscriptRef38,
            final String knownExonUpRangeRef38, final String knownExonDownRangeRef38, final String overridesRef38)
    {
        Type = type;
        FiveGene = fiveGene;
        ThreeGene = threeGene;
        CancerTypes = cancerTypes;
        PubMedId = pubMedId;
        KnownExonTranscript37 = knownExonTranscript37;
        KnownExonUpRange37 = knownExonUpRange37;
        KnownExonDownRange37 = knownExonDownRange37;
        HighImpactPromiscuous = highImpactPromiscuous;
        OverridesRef37 = overridesRef37;
        KnownExonTranscriptRef38 = knownExonTranscriptRef38;
        KnownExonUpRangeRef38 = knownExonUpRangeRef38;
        KnownExonDownRangeRef38 = knownExonDownRangeRef38;
        OverridesRef38 = overridesRef38;
    }

    public boolean isDuplicate(final FusionRefData other)
    {
        if(Type != other.Type)
        {
            return false;
        }

        if(!FiveGene.equals(other.FiveGene) || !ThreeGene.equals(other.ThreeGene))
        {
            return false;
        }

        boolean equalKnownExonRange37 = KnownExonTranscript37.equals(other.KnownExonTranscript37)
                && KnownExonUpRange37.equals(other.KnownExonUpRange37)
                && KnownExonDownRange37.equals(other.KnownExonDownRange37);

        boolean equalKnownExonRange38 = KnownExonTranscriptRef38.equals(other.KnownExonTranscriptRef38)
                && KnownExonUpRangeRef38.equals(other.KnownExonUpRangeRef38)
                && KnownExonDownRangeRef38.equals(other.KnownExonDownRangeRef38);

        return equalKnownExonRange37 || equalKnownExonRange38;
    }

    public String toString()
    {
        return format("%s: %s-%s", Type, FiveGene, ThreeGene);
    }
}
