package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.utils.Nucleotide;

public final class NucleotideFragmentQualEnrichment
{
    public static List<Fragment> qualityFilterFragments(
            int minEvidenceSupport, double minEvidenceFactor, double minHighQualEvidenceFactor, final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        // filter fragments so that each nucleotide has at least 1 base at or above the min-qual threshold, and
        // X fragments (minEvidence) at that base with any qual
        SequenceCount highQualCounts = SequenceCount.nucleotides(minEvidenceSupport, minHighQualEvidenceFactor, highQualFrags);
        SequenceCount rawCounts = SequenceCount.nucleotides(minEvidenceSupport, minEvidenceFactor, fragments);

        return fragments.stream().map(x -> applyQualityFilter(x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private static Fragment applyQualityFilter(final Fragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount)
    {
        // checks whether all nucleotides have qual above the required level - if so return this fragment unch, otherwise build a
        // new fragment just with these filtered loci
        SortedMap<Integer, Nucleotide> nucleotidesByLoci = fragment.nucleotidesByLoci();
        boolean allPresent = true;
        final List<Nucleotide> filteredNucleotides = Lists.newArrayListWithExpectedSize(nucleotidesByLoci.size());
        for(Map.Entry<Integer, Nucleotide> entry : nucleotidesByLoci.entrySet())
        {
            int locus = entry.getKey();
            Nucleotide nucleotide = entry.getValue();
            String fragmentNucleotide = entry.getValue().bases();

            List<String> highQualitySequences = highQualityCount.getMinEvidenceSequences(locus, DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR);
            List<String> rawSequences = rawCount.getMinEvidenceSequences(locus, DEFAULT_MIN_EVIDENCE_FACTOR);
            List<String> allowedSequences = highQualitySequences.stream().filter(x -> rawSequences.contains(x)).collect(Collectors.toList());

            if(allowedSequences.contains(fragmentNucleotide))
            {
                filteredNucleotides.add(nucleotide);
            }
            else
            {
                allPresent = false;
            }
        }

        if(allPresent)
            return fragment;

        Fragment newFragment = new Fragment(
                fragment.reads().get(0), fragment.readGene(), fragment.genes(), filteredNucleotides);

        newFragment.addReads(fragment);

        return newFragment;
    }
}
