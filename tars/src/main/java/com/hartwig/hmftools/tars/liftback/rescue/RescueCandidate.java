package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Input DTO for JunctionRescueResolver, expressed as plain values so the resolver is testable without SAMRecord.
// readBases is only needed for the ref-verify path. mateHintIntrons, when set, biases junction snap toward the
// partner mate's already-rescued junctions.
// NOTE: the generated equals/hashCode compare readBases by array identity, not content; do not use this record
// as a Set/Map key or rely on value equality (it is only ever passed to resolve() as a transient input).
public record RescueCandidate(
        String chromosome, boolean forwardStrand, int readLength, int primaryStart, String primaryCigar,
        int primaryMapq, List<RescueSupplementary> supplementaries, byte[] readBases,
        List<ChrBaseRegion> mateHintIntrons)
{
    public RescueCandidate
    {
        if(mateHintIntrons == null)
            mateHintIntrons = Collections.emptyList();
    }

    public RescueCandidate(
            final String chromosome, final boolean forwardStrand, final int readLength,
            final int primaryStart, final String primaryCigar, final int primaryMapq,
            final List<RescueSupplementary> supplementaries)
    {
        this(chromosome, forwardStrand, readLength, primaryStart, primaryCigar, primaryMapq,
                supplementaries, null, Collections.emptyList());
    }
}
