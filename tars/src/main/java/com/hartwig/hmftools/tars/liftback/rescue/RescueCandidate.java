package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.Collections;
import java.util.List;

// Input DTO for JunctionRescueResolver, expressed as plain values so the resolver is testable without SAMRecord.
// ReadBases is only needed for the ref-verify path. MateHintIntrons, when set, biases junction snap toward the
// partner mate's already-rescued junctions.
public class RescueCandidate
{
    public final String Chromosome;
    public final boolean ForwardStrand;
    public final int ReadLength;
    public final int PrimaryStart;     // 1-based inclusive
    public final String PrimaryCigar;
    public final int PrimaryMapq;
    public final List<RescueSupplementary> Supplementaries;
    public final byte[] ReadBases;     // nullable; only needed for ref-verify path
    public final List<ChrIntron> MateHintIntrons;

    public RescueCandidate(
            final String chromosome, final boolean forwardStrand, final int readLength,
            final int primaryStart, final String primaryCigar, final int primaryMapq,
            final List<RescueSupplementary> supplementaries)
    {
        this(chromosome, forwardStrand, readLength, primaryStart, primaryCigar, primaryMapq,
                supplementaries, null, Collections.emptyList());
    }

    public RescueCandidate(
            final String chromosome, final boolean forwardStrand, final int readLength,
            final int primaryStart, final String primaryCigar, final int primaryMapq,
            final List<RescueSupplementary> supplementaries,
            final byte[] readBases, final List<ChrIntron> mateHintIntrons)
    {
        Chromosome = chromosome;
        ForwardStrand = forwardStrand;
        ReadLength = readLength;
        PrimaryStart = primaryStart;
        PrimaryCigar = primaryCigar;
        PrimaryMapq = primaryMapq;
        Supplementaries = supplementaries;
        ReadBases = readBases;
        MateHintIntrons = mateHintIntrons != null ? mateHintIntrons : Collections.emptyList();
    }
}
