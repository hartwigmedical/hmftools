package com.hartwig.hmftools.lilac.evidence;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment;
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles;
import com.sun.tools.javac.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public final class ExtendEvidence
{
    private final LilacConfig mConfig;
    private final List<Integer> mHeterozygousLoci;
    private final List<AminoAcidFragment> mAminoAcidFragments;
    private final ExpectedAlleles mExpectedAlleles;

    public ExtendEvidence(@NotNull LilacConfig config, @NotNull List<Integer> heterozygousLoci,
            @NotNull List<AminoAcidFragment> aminoAcidFragments, @NotNull ExpectedAlleles expectedAlleles)
    {
        mConfig = config;
        mHeterozygousLoci = heterozygousLoci;
        mAminoAcidFragments = aminoAcidFragments;
        mExpectedAlleles = expectedAlleles;
    }

    public final List<PhasedEvidence> pairedEvidence()
    {
        return Lists.newArrayList();

        /*
        int n;
        int n2;
        List result = new ArrayList();
        if(mConfig.getDebugPhasing())
        {
            PhasedEvidenceFactory.Companion.getLogger().info("    Producing paired evidence");
        }
        if((n2 = 0) <= (n = mHeterozygousLoci.size() - 2))
        {
            void i;
            do
            {
                void $receiver$iv$iv;
                Iterable $receiver$iv;
                List indices = CollectionsKt.listOf((Object[]) new Integer[] { mHeterozygousLoci.get((int) (++i)),
                        mHeterozygousLoci.get((int) (i + true)) });
                Iterable iterable = $receiver$iv = (Iterable) mAminoAcidFragments;
                Collection destination$iv$iv = new ArrayList();
                for(Object element$iv$iv : $receiver$iv$iv)
                {
                    AminoAcidFragment it = (AminoAcidFragment) element$iv$iv;
                    boolean bl = false;
                    if(!it.containsAll(indices))
                    {
                        continue;
                    }
                    destination$iv$iv.add(element$iv$iv);
                }
                List filteredFragments = (List) destination$iv$iv;
                $receiver$iv = filteredFragments;
                if(!(!$receiver$iv.isEmpty()))
                {
                    continue;
                }
                int minTotalFragments = minTotalFragments(indices);
                PhasedEvidence
                        left = PhasedEvidence.Companion.evidence(mAminoAcidFragments, ((Number) indices.get(0)).intValue());
                PhasedEvidence
                        right = PhasedEvidence.Companion.evidence(mAminoAcidFragments, ((Number) indices.get(1)).intValue());
                int[] nArray = CollectionsKt.toIntArray((Collection) indices);
                PhasedEvidence combinedEvidence =
                        PhasedEvidence.Companion.evidence(mAminoAcidFragments, Arrays.copyOf(nArray, nArray.length))
                                .removeSingles(mConfig.getMinFragmentsToRemoveSingle());
                if(combinedEvidence.totalEvidence() >= minTotalFragments
                        && CombineEvidence.INSTANCE.canCombine(left, combinedEvidence, right))
                {
                    result.add(combinedEvidence);
                    if(!mConfig.getDebugPhasing())
                    {
                        continue;
                    }
                    PhasedEvidenceFactory.Companion.getLogger().info("    Paired Evidence: " + combinedEvidence);
                    continue;
                }
                if(!mConfig.getDebugPhasing())
                {
                    continue;
                }
                PhasedEvidenceFactory.Companion.getLogger().info("    FAILED Paired Evidence: " + combinedEvidence);
            } while(i != n);
        }
        return CollectionsKt.sorted((Iterable) result);

         */
    }

    public Pair<PhasedEvidence, Set<PhasedEvidence>> merge(@NotNull PhasedEvidence current, @NotNull Set<PhasedEvidence> others)
    {
        return null;

        /*
        Object left;
        void $receiver$iv$iv;
        void $receiver$iv$iv2;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull((Object) current, (String) "current");
        Intrinsics.checkParameterIsNotNull(others, (String) "others");
        int[] existingIndices = current.getAminoAcidIndices();
        Integer n = ArraysKt.min((int[]) existingIndices);
        if(n == null)
        {
            Intrinsics.throwNpe();
        }
        int minExisting = n;
        Integer n2 = ArraysKt.max((int[]) existingIndices);
        if(n2 == null)
        {
            Intrinsics.throwNpe();
        }
        int maxExisting = n2;
        Iterable iterable = $receiver$iv = (Iterable) others;
        Iterable destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            PhasedEvidence it = (PhasedEvidence) element$iv$iv;
            boolean bl = false;
            if(!(Intrinsics.areEqual((Object) it, (Object) current) ^ true
                    && ArraysKt.contains((int[]) it.getAminoAcidIndices(), (int) maxExisting)))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List othersContainingMax = (List) destination$iv$iv;
        Iterable $receiver$iv2 = others;
        destination$iv$iv = $receiver$iv2;
        Collection destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            PhasedEvidence it = (PhasedEvidence) element$iv$iv;
            boolean bl = false;
            if(!(Intrinsics.areEqual((Object) it, (Object) current) ^ true
                    && ArraysKt.contains((int[]) it.getAminoAcidIndices(), (int) minExisting)))
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        List othersContainingMin = (List) destination$iv$iv2;
        $receiver$iv2 = othersContainingMin;
        if(!$receiver$iv2.isEmpty() && !($receiver$iv2 = (Collection) othersContainingMax).isEmpty())
        {
            Pair<PhasedEvidence, Set<PhasedEvidence>> result;
            left = (PhasedEvidence) othersContainingMin.get(0);
            PhasedEvidence right = (PhasedEvidence) othersContainingMax.get(0);
            return ((PhasedEvidence) left).totalEvidence() > right.totalEvidence()
                    ? (((Set) (result = merge(current, (PhasedEvidence) left, current)).getSecond()).isEmpty()
                    ? merge(current, current, right)
                    : result)
                    : (((Set) (result = merge(current, current, right)).getSecond()).isEmpty()
                            ? merge(current, (PhasedEvidence) left, current)
                            : result);
        }
        left = othersContainingMin;
        if(!left.isEmpty())
        {
            left = (PhasedEvidence) othersContainingMin.get(0);
            return merge(current, (PhasedEvidence) left, current);
        }
        left = othersContainingMax;
        if(!left.isEmpty())
        {
            PhasedEvidence right = (PhasedEvidence) othersContainingMax.get(0);
            return merge(current, current, right);
        }
        return new Pair((Object) current, (Object) SetsKt.emptySet());

         */
    }

    private final int minTotalFragments(List<Integer> indices)
    {
        return mExpectedAlleles.expectedAlleles(indices) * mConfig.MinFragmentsPerAllele;
    }

    private final Pair<PhasedEvidence, Set<PhasedEvidence>> merge(
            PhasedEvidence current, PhasedEvidence left, PhasedEvidence right)
    {
        return null;

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        int[] leftTail = left.unambiguousTailIndices();
        int[] rightHead = right.unambiguousHeadIndices();
        List mergeIndices = CollectionsKt.sorted((Iterable) ArraysKt.distinct((int[]) ArraysKt.plus((int[]) leftTail, (int[]) rightHead)));
        Iterable iterable = $receiver$iv = (Iterable) mAminoAcidFragments;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            AminoAcidFragment it = (AminoAcidFragment) element$iv$iv;
            boolean bl = false;
            if(!it.containsAll(mergeIndices))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List filteredFragments = (List) destination$iv$iv;
        $receiver$iv = filteredFragments;
        if(!$receiver$iv.isEmpty())
        {
            PhasedEvidence combined;
            int minTotalFragments = minTotalFragments(mergeIndices);
            int[] nArray = CollectionsKt.toIntArray((Collection) mergeIndices);
            PhasedEvidence
                    mergeEvidence = PhasedEvidence.Companion.evidence(filteredFragments, Arrays.copyOf(nArray, nArray.length))
                    .removeSingles(mConfig.getMinFragmentsToRemoveSingle());
            if(CombineEvidence.INSTANCE.canCombine(left, mergeEvidence, right)
                    && (combined = CombineEvidence.INSTANCE.combine(left, mergeEvidence, right)).totalEvidence() >= minTotalFragments)
            {
                return new Pair((Object) combined, (Object) SetsKt.setOf((Object[]) new PhasedEvidence[] { left, right }));
            }
        }
        return new Pair((Object) current, (Object) SetsKt.emptySet());

         */
    }

}
