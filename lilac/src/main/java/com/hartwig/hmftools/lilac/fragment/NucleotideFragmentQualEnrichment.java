package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_HIGH_QUAL_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline.RAW_REF_NUCLEOTIDE_COUNTS;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

import org.jetbrains.annotations.Nullable;

public final class NucleotideFragmentQualEnrichment
{
    private NucleotideFragmentQualEnrichment() {}

    public static List<Fragment> qualityFilterFragments(
            final HlaContext context, final Collection<Fragment> fragments, final Iterable<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        // filter fragments so that each nucleotide has at least a VAF above some threshold at or above the min-qual threshold, and
        // VAF above some threshold (minEvidenceFactor) at that base with any qual, also automatically include nucleotidesthat have low
        // raw depth
        SequenceCount highQualCounts = SequenceCount.buildFromNucleotides(MIN_HIGH_QUAL_EVIDENCE_FACTOR, highQualFrags);
        SequenceCount rawCounts = SequenceCount.buildFromNucleotides(MIN_EVIDENCE_FACTOR, fragments);

        SequenceCount rawNucleotideCounts = RAW_REF_NUCLEOTIDE_COUNTS.get(context.Gene);
        return fragments.stream()
                .map(x -> applyQualityFilter(x, highQualCounts, rawCounts, rawNucleotideCounts))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    public static Fragment applyQualityFilter(final Fragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount,
            @Nullable final SequenceCount rawNucleotideCounts)
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
            if(rawNucleotideCounts != null && rawNucleotideCounts.get(locus).size() < MIN_DEPTH_FILTER)
            {
                filteredNucleotides.add(nucleotide);
                continue;
            }

            String fragmentNucleotide = entry.getValue().bases();

            Set<String> highQualitySequences = Sets.newHashSet(highQualityCount.getMinEvidenceSequences(locus));
            Set<String> rawSequences = Sets.newHashSet(rawCount.getMinEvidenceSequences(locus));
            Set<String> allowedSequences = Sets.newHashSet(highQualitySequences);
            allowedSequences.retainAll(rawSequences);
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
