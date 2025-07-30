package com.hartwig.hmftools.lilac.fragment;

import java.util.List;
import java.util.NavigableSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
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

        NavigableSet<Integer> loci = Sets.newTreeSet();
        loci.addAll(referenceCount.seqCountsByLoci().keySet());
        loci.addAll(tumorCount.seqCountsByLoci().keySet());

        for(int locus : loci)
        {
            boolean inRef = referenceCount.seqCountsByLoci().containsKey(locus);
            boolean inTum = tumorCount.seqCountsByLoci().containsKey(locus);
            Multiset<String> refSeqCounts = referenceCount.get(locus);
            Multiset<String> tumorSeqCounts = tumorCount.get(locus);
            if(inRef)
            {
                int refDepth = referenceCount.depth(locus);
                final int locusIndex = locus;
                refSeqCounts.entrySet().stream()
                        .filter(x -> !tumorSeqCounts.contains(x.getElement()))
                        .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x.getElement(), x.getCount(), refDepth, 0, 0)));
            }

            if(inTum)
            {
                int tumorDepth = tumorCount.depth(locus);
                final int locusIndex = locus;
                tumorSeqCounts.entrySet().stream()
                        .filter(x -> !refSeqCounts.contains(x.getElement()))
                        .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(locusIndex, x.getElement(), 0, 0, x.getCount(), tumorDepth)));
            }
        }

        return seqCountDiffs;
    }
}
