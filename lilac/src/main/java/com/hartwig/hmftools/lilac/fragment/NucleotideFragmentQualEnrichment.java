package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;

public final class NucleotideFragmentQualEnrichment
{
    public static List<Fragment> qualityFilterFragments(final HlaContext context, int minEvidenceSupport, double minEvidenceFactor,
            double minHighQualEvidenceFactor, final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        // filter fragments so that each nucleotide has at least a VAF above some threshold at or above the min-qual threshold, and
        // VAF above some threshold (minEvidenceFactor) at that base with any qual, also automatically include nucleotidesthat have low
        // raw depth
        SequenceCount highQualCounts = SequenceCount.nucleotides(minEvidenceSupport, minHighQualEvidenceFactor, highQualFrags);
        SequenceCount rawCounts = SequenceCount.nucleotides(minEvidenceSupport, minEvidenceFactor, fragments);

        return fragments.stream().map(x -> applyQualityFilter(context, x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private static Fragment applyQualityFilter(final HlaContext context, final Fragment fragment, final SequenceCount highQualityCount,
            final SequenceCount rawCount)
    {
        // checks whether all nucleotides have qual above the required level - if so return this fragment, otherwise build a
        // new fragment just with these filtered loci
        SortedMap<Integer, Nucleotide> nucleotidesByLoci = fragment.nucleotidesByLoci();
        boolean allPresent = true;
        final List<Nucleotide> filteredNucleotides = Lists.newArrayListWithExpectedSize(nucleotidesByLoci.size());
        for(Map.Entry<Integer, Nucleotide> entry : nucleotidesByLoci.entrySet())
        {
            int locus = entry.getKey();
            Nucleotide nucleotide = entry.getValue();
            String fragmentNucleotide = entry.getValue().bases();

            Set<String> highQualitySequences = Sets.newHashSet(
                    highQualityCount.getMinEvidenceSequences(locus, DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR));
            Set<String> rawSequences = Sets.newHashSet(rawCount.getMinEvidenceSequences(locus, DEFAULT_MIN_EVIDENCE_FACTOR));
            Set<String> lowDepthSequences = rawCount.getLowRawDepthSequences(context.geneName(), locus, DEFAULT_MIN_DEPTH_FILTER);

            Set<String> allowedSequences = Sets.newHashSet(highQualitySequences);
            allowedSequences.retainAll(rawSequences);
            allowedSequences.addAll(lowDepthSequences);

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
