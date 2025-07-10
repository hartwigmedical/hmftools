package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Integer.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class SequenceCountDiff
{
    public final int Loci;
    public final String Sequence;
    public final int ReferenceCount;
    public final int ReferenceDepth;
    public final int TumorCount;
    public final int TumorDepth;

    public SequenceCountDiff(int locus, final String sequence, int referenceCount, int referenceDepth, int tumorCount, int tumorDepth)
    {
        Loci = locus;
        Sequence = sequence;
        ReferenceCount = referenceCount;
        ReferenceDepth = referenceDepth;
        TumorCount = tumorCount;
        TumorDepth = tumorDepth;
    }

    public static List<SequenceCountDiff> create(final SequenceCount referenceCount, final SequenceCount tumorCount)
    {
        final List<SequenceCountDiff> seqCountDiffs = Lists.newArrayList();

        int minLength = min(referenceCount.getLength(), tumorCount.getLength());

        for(int locus = 0; locus < minLength; ++locus) 
        {
            int refDepth = referenceCount.depth(locus);
            int tumorDepth = tumorCount.depth(locus);

            Map<String, Integer> refSeqCounts = referenceCount.get(locus);
            Map<String, Integer> tumorSeqCounts = tumorCount.get(locus);

            final int locusIndex = locus;

            refSeqCounts.keySet().stream()
                    .filter(x -> !tumorSeqCounts.containsKey(x))
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x, refSeqCounts.get(x), refDepth, 0, tumorDepth)));

            tumorSeqCounts.keySet().stream()
                    .filter(x -> !refSeqCounts.containsKey(x))
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x, 0, refDepth, tumorSeqCounts.get(x), tumorDepth)));
        }

        for(int locus = minLength; locus < tumorCount.getLength(); ++locus) 
        {
            int tumorDepth = tumorCount.depth(locus);

            Map<String, Integer> tumorSeqCounts = tumorCount.get(locus);

            final int locusIndex = locus;

            tumorSeqCounts.keySet().stream()
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x, 0, 0, tumorSeqCounts.get(x), tumorDepth)));
        }

        for(int locus = minLength; locus < referenceCount.getLength(); ++locus) 
        {
            int refDepth = referenceCount.depth(locus);

            Map<String, Integer> refSeqCounts = referenceCount.get(locus);

            final int locusIndex = locus;

            refSeqCounts.keySet().stream()
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x, refSeqCounts.get(x), refDepth, 0, 0)));
        }

        return seqCountDiffs;
    }
}
