package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class NucleotideSpliceEnrichment
{
    private final Set<Integer> mAminoAcidBoundary;

    public NucleotideSpliceEnrichment(final Set<Integer> aminoAcidBoundary)
    {
        mAminoAcidBoundary = aminoAcidBoundary;
    }

    public List<Fragment> applySpliceInfo(final Iterable<Fragment> fragments, final Iterable<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        SequenceCount nucleotideCounts = SequenceCount.buildFromNucleotides(MIN_EVIDENCE_FACTOR, highQualFrags);
        Set<Integer> nucleotideExonBoundaryStarts = mAminoAcidBoundary.stream().map(x -> x * 3).collect(Collectors.toSet());
        List<Integer> homLoci = Lists.newArrayList(nucleotideCounts.homozygousLoci());

        Set<Integer> homStarts = nucleotideExonBoundaryStarts.stream().filter(homLoci::contains).collect(Collectors.toSet());

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
                    addEnd(fragment, homEnd, nucleotideCounts);
            }

            results.add(fragment);
        }

        return results;
    }

    private static boolean missingStart(int index, final Fragment fragment)
    {
        return !fragment.containsNucleotideLocus(index)
                && fragment.containsNucleotideLocus(index + 1)
                && fragment.containsNucleotideLocus(index + 2);
    }

    private static boolean missingEnd(int index, final Fragment fragment)
    {
        return fragment.containsNucleotideLocus(index)
                && !fragment.containsNucleotideLocus(index + 1)
                && !fragment.containsNucleotideLocus(index + 2);
    }

    private static void addStart(final Fragment fragment, int index, final SequenceCount nucleotideCounts)
    {
        fragment.addHighQualNucleotide(index, nucleotideCounts.getMinEvidenceSequences(index).get(0));
    }

    private static void addEnd(final Fragment fragment, int index, final SequenceCount nucleotideCounts)
    {
        fragment.addHighQualNucleotide(index + 1, nucleotideCounts.getMinEvidenceSequences(index + 1).get(0));
        fragment.addHighQualNucleotide(index + 2, nucleotideCounts.getMinEvidenceSequences(index + 2).get(0));
    }
}
