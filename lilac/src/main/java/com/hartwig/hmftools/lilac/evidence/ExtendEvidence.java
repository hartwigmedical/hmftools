package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.listMax;
import static com.hartwig.hmftools.lilac.LilacUtils.listMin;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import org.apache.commons.math3.util.Pair;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public final class ExtendEvidence
{
    private final LilacConfig mConfig;
    private final List<Integer> mHeterozygousLoci;
    private final List<AminoAcidFragment> mAminoAcidFragments;
    private final ExpectedAlleles mExpectedAlleles;

    public ExtendEvidence(final LilacConfig config, final List<Integer> heterozygousLoci,
            final List<AminoAcidFragment> aminoAcidFragments, final ExpectedAlleles expectedAlleles)
    {
        mConfig = config;
        mHeterozygousLoci = heterozygousLoci;
        mAminoAcidFragments = aminoAcidFragments;
        mExpectedAlleles = expectedAlleles;
    }

    public final List<PhasedEvidence> pairedEvidence()
    {
        List<PhasedEvidence> results = Lists.newArrayList();

        if (mConfig.DebugPhasing)
        {
            LL_LOGGER.info("  producing paired evidence");
        }

        for(int i = 0; i < mHeterozygousLoci.size() - 1; ++i)
        {
            List<Integer> indices = Lists.newArrayList(mHeterozygousLoci.get(i), mHeterozygousLoci.get(i + 1));

            final List<AminoAcidFragment> filteredFragments = mAminoAcidFragments.stream()
                    .filter(x -> x.containsAll(indices)).collect(Collectors.toList());

            if(!filteredFragments.isEmpty())
            {
                int minTotalFragments = minTotalFragments(indices);

                PhasedEvidence left = PhasedEvidence.evidence(mAminoAcidFragments, Lists.newArrayList(indices.get(0)));
                PhasedEvidence right = PhasedEvidence.evidence(mAminoAcidFragments, Lists.newArrayList(indices.get(1)));

                PhasedEvidence combinedEvidence = PhasedEvidence.evidence(mAminoAcidFragments, indices)
                        .removeSingles(mConfig.MinFragmentsToRemoveSingle);

                if (combinedEvidence.totalEvidence() >= minTotalFragments && CombineEvidence.canCombine(left, combinedEvidence, right))
                {
                    results.add(combinedEvidence);

                    if (mConfig.DebugPhasing)
                    {
                        LL_LOGGER.info("  paired Evidence: {}", combinedEvidence);
                    }
                }
                else
                {
                    if (mConfig.DebugPhasing)
                    {
                        LL_LOGGER.info("  failed paired Evidence: {}", combinedEvidence);
                    }
                }
            }
        }

        Collections.sort(results);
        return results;
    }

    public Pair<PhasedEvidence, Set<PhasedEvidence>> merge(final PhasedEvidence current, final Set<PhasedEvidence> others)
    {
        // List<Integer> existingIndices = current.getAminoAcidIndexList();

        int minExisting = listMin(current.getAminoAcidIndexList());
        int maxExisting = listMax(current.getAminoAcidIndexList());

        List<PhasedEvidence> othersContainingMax = others.stream()
                .filter(x -> x != current && x.getAminoAcidIndexList().contains(maxExisting)).collect(Collectors.toList());

        List<PhasedEvidence> othersContainingMin = others.stream()
                .filter(x -> x != current && x.getAminoAcidIndexList().contains(minExisting)).collect(Collectors.toList());

        if (!othersContainingMin.isEmpty() && !othersContainingMax.isEmpty())
        {
            PhasedEvidence left = othersContainingMin.get(0);
            PhasedEvidence right = othersContainingMax.get(0);

            if (left.totalEvidence() > right.totalEvidence())
            {
                Pair<PhasedEvidence, Set<PhasedEvidence>> result = merge(current, left, current);

                if (result.getSecond().isEmpty())
                    return merge(current, current, right);
                else
                    return result;
            }
            else
            {
                Pair<PhasedEvidence, Set<PhasedEvidence>> result = merge(current, current, right);

                if (result.getSecond().isEmpty())
                    return merge(current, left, current);
                else
                    return result;
            }
        }

        if (!othersContainingMin.isEmpty())
        {
            PhasedEvidence left = othersContainingMin.get(0);
            return merge(current, left, current);
        }

        if (!othersContainingMax.isEmpty())
        {
            PhasedEvidence right = othersContainingMax.get(0);
            return merge(current, current, right);
        }

        return new Pair(current, Sets.newHashSet());
    }

    private final int minTotalFragments(List<Integer> indices)
    {
        return mExpectedAlleles.expectedAlleles(indices) * mConfig.MinFragmentsPerAllele;
    }

    private final Pair<PhasedEvidence, Set<PhasedEvidence>> merge( PhasedEvidence current, PhasedEvidence left, PhasedEvidence right)
    {
        List<Integer> leftTail = left.unambiguousTailIndices();
        List<Integer> rightHead = right.unambiguousHeadIndices();
        List<Integer> mergeIndices = leftTail.stream().collect(Collectors.toList());

        rightHead.stream().filter(x -> !leftTail.contains(x)).forEach(x -> mergeIndices.add(x));
        Collections.sort(mergeIndices);

        List<AminoAcidFragment> filteredFragments = mAminoAcidFragments.stream()
                .filter(x -> x.containsAll(mergeIndices)).collect(Collectors.toList());

        if (!filteredFragments.isEmpty())
        {
            int minTotalFragments = minTotalFragments(mergeIndices);

            PhasedEvidence mergeEvidence = PhasedEvidence.evidence(filteredFragments, mergeIndices).removeSingles(mConfig.MinFragmentsToRemoveSingle);

            if (CombineEvidence.canCombine(left, mergeEvidence, right))
            {
                PhasedEvidence combined = CombineEvidence.combine(left, mergeEvidence, right);

                if (combined.totalEvidence() >= minTotalFragments)
                {
                    return new Pair(combined, Sets.newHashSet(left, right));
                }
            }
        }

        return new Pair(current, Sets.newHashSet());
    }
}
