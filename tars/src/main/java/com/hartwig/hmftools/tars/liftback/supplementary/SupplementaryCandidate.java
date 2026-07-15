package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Input DTO for SupplementaryResolver, expressed as plain values so the resolver is testable without SAMRecord.
// readBases is only needed for the ref-verify path. mateHintIntrons, when set, biases junction snap toward the
// partner mate's already-resolved junctions.
// NOTE: the generated equals/hashCode compare readBases by array identity, not content; do not use this record
// as a Set/Map key or rely on value equality (it is only ever passed to resolve() as a transient input).
public record SupplementaryCandidate(
        String chromosome, boolean forwardStrand, int readLength, int primaryStart, String primaryCigar,
        int primaryMapq, List<SupplementaryRecord> supplementaries, byte[] readBases,
        List<ChrBaseRegion> mateHintIntrons)
{
    public SupplementaryCandidate
    {
        if(mateHintIntrons == null)
        {
            mateHintIntrons = Collections.emptyList();
        }
    }

    public SupplementaryCandidate(
            final String chromosome, final boolean forwardStrand, final int readLength,
            final int primaryStart, final String primaryCigar, final int primaryMapq,
            final List<SupplementaryRecord> supplementaries)
    {
        this(chromosome, forwardStrand, readLength, primaryStart, primaryCigar, primaryMapq,
                supplementaries, null, Collections.emptyList());
    }
}
