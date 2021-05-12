package com.hartwig.hmftools.lilac.candidates;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class AminoAcidFiltering
{
    private final List<Integer> mAminoAcidBoundaries;

    public AminoAcidFiltering(final List<Integer> aminoAcidBoundaries)
    {
        mAminoAcidBoundaries = aminoAcidBoundaries;
    }

    public List<HlaSequenceLoci> aminoAcidCandidates(final List<HlaSequenceLoci> candidates, final SequenceCount aminoAcidCount)
    {
        List<HlaSequenceLoci> results = Lists.newArrayList();
        results.addAll(candidates);

        Set<Integer> locations = Sets.newHashSet();

        for(int loci = 0; loci < aminoAcidCount.getLength(); ++loci)
        {
            if(mAminoAcidBoundaries.contains(loci))
                continue;

            // int depth = aminoAcidCount.depth(loci);

            List<String> expectedSequences = aminoAcidCount.getMinCountSequences(loci);

            final int lociConst = loci;

            results = results.stream()
                    .filter(x -> x.consistentWithAny(expectedSequences, Lists.newArrayList(lociConst)))
                    .collect(Collectors.toList());
        }

        return results;
    }

}
