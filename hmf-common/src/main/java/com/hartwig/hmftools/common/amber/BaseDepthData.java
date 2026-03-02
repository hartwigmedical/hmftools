package com.hartwig.hmftools.common.amber;

public record BaseDepthData(AmberBase ref,
                            AmberBase alt,
                            int readDepth,
                            int indelCount,
                            int refSupport,
                            int altSupport)
{
}
