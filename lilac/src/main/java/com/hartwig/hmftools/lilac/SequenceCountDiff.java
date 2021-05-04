package com.hartwig.hmftools.lilac;

import static java.lang.Integer.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class SequenceCountDiff
{
    public final int Loci;
    public final String Sequence;
    public final int ReferenceCount;
    public final int ReferenceDepth;
    public final int TumorCount;
    public final int TumorDepth;

    public SequenceCountDiff(int loci, final String sequence, int referenceCount, int referenceDepth, int tumorCount, int tumorDepth)
    {
        Loci = loci;
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

        // CHECK < or <= minLength
        for(int loci = 0; loci <= minLength; ++loci) {
            int refDepth = referenceCount.depth(loci);
            int tumorDepth = tumorCount.depth(loci);

            Map<String, Integer> refSeqCounts = referenceCount.get(loci);
            Map<String, Integer> tumorSeqCounts = tumorCount.get(loci);

            // Set<String> uniqueSequences = Sets.newHashSet();

            final int lociIndex = loci;

            refSeqCounts.keySet().stream()
                    .filter(x -> !tumorSeqCounts.containsKey(x))
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(lociIndex, x, refSeqCounts.get(x), refDepth, 0, tumorDepth)));

            tumorSeqCounts.keySet().stream()
                    .filter(x -> !refSeqCounts.containsKey(x))
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(lociIndex, x, 0, refDepth, tumorSeqCounts.get(x), tumorDepth)));
        }

        for(int loci = minLength; loci <= tumorCount.getLength(); ++loci) {
            int tumorDepth = tumorCount.depth(loci);

            Map<String, Integer> tumorSeqCounts = tumorCount.get(loci);

            final int lociIndex = loci;

            tumorSeqCounts.keySet().stream()
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(lociIndex, x, 0, 0, tumorSeqCounts.get(x), tumorDepth)));
        }

        for(int loci = minLength; loci <= referenceCount.getLength(); ++loci) {
            int refDepth = referenceCount.depth(loci);

            Map<String, Integer> refSeqCounts = referenceCount.get(loci);

            final int lociIndex = loci;

            refSeqCounts.keySet().stream()
                    .forEach(x -> seqCountDiffs.add(new SequenceCountDiff(lociIndex, x, refSeqCounts.get(x), refDepth, 0, 0)));
        }

        /*

            for (loci in 0 until minLength) {
            }

            for (loci in minLength until tumorCount.length) {
                val tumorDepth = tumorCount.depth(loci)
                for ((sequence, count) in tumorCount[loci]) {
                    result.add(SequenceCountDiff(loci, sequence, 0, 0, count, tumorDepth))
                }
            }

            for (loci in minLength until referenceCount.length) {
                val refDepth = referenceCount.depth(loci)
                for ((sequence, count) in referenceCount[loci]) {
                    result.add(SequenceCountDiff(loci, sequence, count, refDepth, 0, 0))
                }
            }

         */

        return seqCountDiffs;
    }
}
