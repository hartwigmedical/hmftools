package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class NucleotideSpliceEnrichment
{
    private final byte mMinBaseQuality;
    private final double mMinEvidence;
    private final Set<Integer> mAminoAcidBoundary;

    public NucleotideSpliceEnrichment(double minEvidence, final Set<Integer> aminoAcidBoundary)
    {
        mMinBaseQuality = LOW_BASE_QUAL_THRESHOLD;
        mMinEvidence = minEvidence;
        mAminoAcidBoundary = aminoAcidBoundary;
    }

    public List<Fragment> applySpliceInfo(final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        SequenceCount nucleotideCounts = SequenceCount.nucleotides(mMinEvidence, highQualFrags);
        Set<Integer> nucleotideExonBoundaryStarts = mAminoAcidBoundary.stream().map(x -> x * 3).collect(Collectors.toSet());
        List<Integer> homLoci = nucleotideCounts.homozygousIndices();

        Set<Integer> homStarts = nucleotideExonBoundaryStarts.stream().filter(x -> homLoci.contains(x)).collect(Collectors.toSet());

        Set<Integer> homEnds = nucleotideExonBoundaryStarts.stream()
                .filter(x -> homLoci.contains(x + 1) && homLoci.contains(x + 2)).collect(Collectors.toSet());

        final List<Fragment> results = Lists.newArrayList();

        for (Fragment fragment : fragments)
        {
            for (Integer homStart : homStarts)
            {
                if (missingStart(homStart, fragment))
                {
                    addStart(fragment, homStart, nucleotideCounts);
                }
            }

            for (Integer homEnd : homEnds)
            {
                if (missingEnd(homEnd, fragment))
                    addEnd(fragment, homEnd, nucleotideCounts);
            }

            results.add(fragment);
        }

        return results;
    }

    private boolean missingStart(int index, Fragment fragment)
    {
        return !fragment.containsNucleotideLocus(index) && fragment.containsAllNucleotideLoci(Lists.newArrayList(index + 1, index + 2));
    }

    private boolean missingEnd(int index, Fragment fragment)
    {
        return fragment.containsNucleotideLocus(index) && !fragment.containsNucleotideLocus(index + 1) && !fragment.containsNucleotideLocus(index + 2);
    }

    private void addStart(final Fragment fragment, int index, SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index, nucleotideCounts.getMinCountSequences(index).get(0), mMinBaseQuality);
    }

    private void addEnd(final Fragment fragment, int index, SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index + 1, nucleotideCounts.getMinCountSequences(index + 1).get(0), mMinBaseQuality);
        fragment.addNucleotide(index + 2, nucleotideCounts.getMinCountSequences(index + 2).get(0), mMinBaseQuality);
    }
}
