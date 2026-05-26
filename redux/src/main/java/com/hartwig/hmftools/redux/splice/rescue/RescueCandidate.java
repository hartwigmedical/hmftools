package com.hartwig.hmftools.redux.splice.rescue;

import java.util.Collections;
import java.util.List;

// Input DTO for JunctionRescueResolver. Carries the per-mate state needed to decide a merge,
// expressed as plain values so the resolver can be tested without involving SAMRecord. Construction
// from real SAM records happens in the adapter that wires the resolver into LiftBackResolver.
//
// ReadBases is optional — only required for the ref-verify path (no-supp rescue via reference
// lookup). MateHintIntrons is optional — when populated with the partner mate's previously-
// rescued junctions, the snap loop will prefer one of those as a tie-breaker.
public class RescueCandidate
{
    public final String Chromosome;
    public final boolean ForwardStrand;
    public final int ReadLength;
    public final int PrimaryStart;     // 1-based
    public final String PrimaryCigar;
    public final int PrimaryMapq;
    public final List<RescueSupplementary> Supplementaries;
    public final byte[] ReadBases;     // nullable
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
