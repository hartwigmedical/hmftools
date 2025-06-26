package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class NucleotideSpliceEnrichment
{
    private final byte mMinBaseQuality;
    private final int mMinValFilterDepth;
    private final double mMinEvidenceFactor;
    private final Set<Integer> mAminoAcidBoundary;

    public NucleotideSpliceEnrichment(final int minVafFilterDepth, final double minEvidenceFactor,
            final Set<Integer> aminoAcidBoundary)
    {
        mMinBaseQuality = LOW_BASE_QUAL_THRESHOLD;
        mMinValFilterDepth = minVafFilterDepth;
        mMinEvidenceFactor = minEvidenceFactor;
        mAminoAcidBoundary = aminoAcidBoundary;
    }

    public List<Fragment> applySpliceInfo(final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        SequenceCount nucleotideCounts = SequenceCount
                .nucleotides(mMinValFilterDepth, mMinEvidenceFactor, highQualFrags);
        Set<Integer> nucleotideExonBoundaryStarts = mAminoAcidBoundary.stream().map(x -> x * 3).collect(Collectors.toSet());
        NavigableSet<Integer> homLoci = nucleotideCounts.homozygousLoci();

        Set<Integer> homStarts = nucleotideExonBoundaryStarts.stream().filter(x -> homLoci.contains(x)).collect(Collectors.toSet());

        Set<Integer> homEnds = nucleotideExonBoundaryStarts.stream()
                .filter(x -> homLoci.contains(x + 1) && homLoci.contains(x + 2)).collect(Collectors.toSet());

        final List<Fragment> results = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            for(Integer homStart : homStarts)
            {
                if(missingStart(homStart, fragment))
                {
                    addStart(fragment, homStart, nucleotideCounts);
                }
            }

            for(Integer homEnd : homEnds)
            {
                if(missingEnd(homEnd, fragment))
                {
                    addEnd(fragment, homEnd, nucleotideCounts);
                }
            }

            results.add(fragment);
        }

        return results;
    }

    private boolean missingStart(final int index, final Fragment fragment)
    {
        return !fragment.containsNucleotideLocus(index) && fragment.containsAllNucleotideLoci(Lists.newArrayList(index + 1, index + 2));
    }

    private boolean missingEnd(final int index, final Fragment fragment)
    {
        return fragment.containsNucleotideLocus(index) && !fragment.containsNucleotideLocus(index + 1) && !fragment.containsNucleotideLocus(
                index + 2);
    }

    private void addStart(final Fragment fragment, final int index, final SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index, nucleotideCounts.getMinEvidenceSequences(index).get(0), mMinBaseQuality);
    }

    private void addEnd(final Fragment fragment, final int index, final SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index + 1, nucleotideCounts.getMinEvidenceSequences(index + 1).get(0), mMinBaseQuality);
        fragment.addNucleotide(index + 2, nucleotideCounts.getMinEvidenceSequences(index + 2).get(0), mMinBaseQuality);
    }
}
