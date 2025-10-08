package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class NucleotideSpliceEnrichment
{
    private final byte mMinBaseQuality;
    private final Optional<Set<Integer>> mAminoAcidBoundary_;

        public NucleotideSpliceEnrichment(final Optional<Set<Integer>> aminoAcidBoundary)
    {
        mMinBaseQuality = LOW_BASE_QUAL_THRESHOLD;
        mAminoAcidBoundary_ = aminoAcidBoundary;
    }

    public List<Fragment> applySpliceInfo(final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        SequenceCount nucleotideCounts = SequenceCount.buildFromNucleotides(MIN_EVIDENCE_FACTOR, highQualFrags);
        Optional<Set<Integer>> nucleotideExonBoundaryStarts_ = Optional.empty();
        if(mAminoAcidBoundary_.isPresent())
        {
            nucleotideExonBoundaryStarts_ = Optional.of(mAminoAcidBoundary_.get().stream().map(x -> x * 3).collect(Collectors.toSet()));
        }

        List<Integer> homLoci = Lists.newArrayList(nucleotideCounts.homozygousLoci());

        // TODO: All of this is only up to max common amino acid boundary.
        Optional<Set<Integer>> homStarts_ = Optional.empty();
        Optional<Set<Integer>> homEnds_ = Optional.empty();
        if(nucleotideExonBoundaryStarts_.isPresent())
        {
            homStarts_ = Optional.of(nucleotideExonBoundaryStarts_.get().stream().filter(x -> homLoci.contains(x)).collect(Collectors.toSet()));
            homEnds_ = Optional.of(nucleotideExonBoundaryStarts_.get().stream()
                    .filter(x -> homLoci.contains(x + 1) && homLoci.contains(x + 2)).collect(Collectors.toSet()));
        }

        final List<Fragment> results = Lists.newArrayList();

        // TODO: this adds boundary start/end based on what we have seen elsewhere.
        for(Fragment fragment : fragments)
        {
            for(Integer homStart : homStarts_.orElse(Sets.newTreeSet()))
            {
                if(missingStart(homStart, fragment))
                {
                    addStart(fragment, homStart, nucleotideCounts);
                }
            }

            for(Integer homEnd : homEnds_.orElse(Sets.newTreeSet()))
            {
                if(missingEnd(homEnd, fragment))
                    addEnd(fragment, homEnd, nucleotideCounts);
            }

            results.add(fragment);
        }

        return results;
    }

    private boolean missingStart(int index, final Fragment fragment)
    {
        return !fragment.containsNucleotideLocus(index)
                && fragment.containsNucleotideLocus(index + 1)
                && fragment.containsNucleotideLocus(index + 2);
    }

    private boolean missingEnd(int index, final Fragment fragment)
    {
        return fragment.containsNucleotideLocus(index)
                && !fragment.containsNucleotideLocus(index + 1)
                && !fragment.containsNucleotideLocus(index + 2);
    }

    private void addStart(final Fragment fragment, int index, final SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index, nucleotideCounts.getMinEvidenceSequences(index).get(0), mMinBaseQuality);
    }

    private void addEnd(final Fragment fragment, int index, final SequenceCount nucleotideCounts)
    {
        fragment.addNucleotide(index + 1, nucleotideCounts.getMinEvidenceSequences(index + 1).get(0), mMinBaseQuality);
        fragment.addNucleotide(index + 2, nucleotideCounts.getMinEvidenceSequences(index + 2).get(0), mMinBaseQuality);
    }
}
