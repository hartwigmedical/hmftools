package com.hartwig.hmftools.common.purple;

import org.jetbrains.annotations.NotNull;

public enum GermlineStatus
{
    HOM_DELETION,
    HET_DELETION,
    AMPLIFICATION,
    NOISE,
    DIPLOID,
    UNKNOWN;

    public static GermlineStatus fromString(@NotNull final String status)
    {
        if(status.equals("SOMATIC"))
        {
            return DIPLOID;
        }

        if(status.equals("CLUSTER"))
        {
            return UNKNOWN;
        }

        if(status.equals("GERMLINE_HOM_DELETION"))
        {
            return HOM_DELETION;
        }

        if(status.equals("GERMLINE_HET_DELETION"))
        {
            return HET_DELETION;
        }

        if(status.equals("GERMLINE_AMPLIFICATION"))
        {
            return AMPLIFICATION;
        }

        if(status.equals("GERMLINE_NOISE"))
        {
            return NOISE;
        }

        return GermlineStatus.valueOf(status);
    }
}
