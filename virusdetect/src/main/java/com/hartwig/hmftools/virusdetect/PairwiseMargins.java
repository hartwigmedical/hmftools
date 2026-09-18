package com.hartwig.hmftools.virusdetect;

import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;

// Per ordered within-group contig pair, how decisively the subject fits shared reads better than the opponent.
// Quantifies the distribution of positive winning alignment margins and the read count.
// This is a key input to the representative contig selection scheme, used to determine if there is significant support
// for a rival virus strain.
public class PairwiseMargins
{
    // Subject -> opponent -> (winning margin in bases -> read count)
    private final Map<ContigPair, NavigableMap<Integer, Integer>> mMarginCounts;
    private final Map<ContigPair, Integer> mSharedReads;
    private final double mMeanReadLength;

    public PairwiseMargins(
            Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts, Map<ContigPair, Integer> sharedReads, double meanReadLength)
    {
        mMarginCounts = marginCounts;
        mSharedReads = sharedReads;
        mMeanReadLength = meanReadLength;
    }

    // Reads whose best alignment fits the subject at least minMargin divergent bases better than the opponent.
    public int challengeReads(ViralContig subject, ViralContig opponent, int minMargin)
    {
        NavigableMap<Integer, Integer> margins = mMarginCounts.get(new ContigPair(subject, opponent));
        if(margins == null)
        {
            return 0;
        }
        return margins.tailMap(minMargin, true).values().stream().mapToInt(Integer::intValue).sum();
    }

    // Reads aligned to both contigs (regardless of which one they favour).
    public int sharedReads(ViralContig subject, ViralContig opponent)
    {
        return mSharedReads.getOrDefault(new ContigPair(subject, opponent), 0);
    }

    public Set<ContigPair> pairs()
    {
        return mSharedReads.keySet();
    }

    public double meanReadLength()
    {
        return mMeanReadLength;
    }

    public record ContigPair(
            ViralContig subject,
            ViralContig opponent
    )
    {
    }
}
