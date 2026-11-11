package com.hartwig.hmftools.viridian.integration;

import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSv;

// A viral integration into the tumor host genome: an SV whose insert aligns to a virus genome.
public record Integration(
    CandidateSv sv,
    ViralInsertAlignment alignment
)
{
}
