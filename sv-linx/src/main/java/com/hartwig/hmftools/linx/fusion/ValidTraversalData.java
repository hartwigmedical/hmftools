package com.hartwig.hmftools.linx.fusion;

import com.hartwig.hmftools.linx.types.LinkedPair;

public class ValidTraversalData
{
    public final LinkedPair Link;
    public final boolean IsValid;
    public final int FusionDirection;
    public final boolean IsPrecodingUpstream;

    public ValidTraversalData(final LinkedPair link, boolean isValid, int fusionDirection, boolean isPrecodingUpstream)
    {
        Link = link;
        IsValid = isValid;
        FusionDirection = fusionDirection;
        IsPrecodingUpstream = isPrecodingUpstream;
    }

    public boolean matches(final LinkedPair link, int fusionDirection, boolean isPrecodingUpstream)
    {
        return (fusionDirection == FusionDirection && isPrecodingUpstream == IsPrecodingUpstream && Link.matches(link));
    }
}
