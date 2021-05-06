package com.hartwig.hmftools.lilac.fragment;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.SequenceCount;

public class AminoAcidQualEnrichment
{
    private final int mMinEvidence;

    public final List<AminoAcidFragment> enrich(final List<NucleotideFragment> nucleotideFragments)
    {
        // TODO
        return Lists.newArrayList();
        /*
        void $receiver$iv$iv;
        NucleotideFragment it;
        Iterable $receiver$iv$iv2;
        Iterable $receiver$iv;
        AminoAcidFragment aminoAcidFragment;
        Collection collection;
        void $receiver$iv$iv3;
        Iterable $receiver$iv2;
        Intrinsics.checkParameterIsNotNull(nucleotideFragments, (String) "nucleotideFragments");
        Iterable iterable = $receiver$iv2 = (Iterable) nucleotideFragments;
        Iterable destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv3)
        {
            Iterator it2;
            NucleotideFragment nucleotideFragment = (NucleotideFragment) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            aminoAcidFragment = ((NucleotideFragment) ((Object) it2)).toAminoAcidFragment();
            collection.add(aminoAcidFragment);
        }
        List qualityFilteredAminoAcidFragments = (List) destination$iv$iv;
        SequenceCount highQualityAminoAcidCounts = SequenceCount.Companion.aminoAcids(this.mMinEvidence, qualityFilteredAminoAcidFragments);
        destination$iv$iv = nucleotideFragments;
        Iterator iterator = $receiver$iv;
        Iterable destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            it = (NucleotideFragment) element$iv$iv;
            boolean bl = false;
            if(!it.isNotEmpty())
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv2;
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv2)
        {
            it = (NucleotideFragment) item$iv$iv;
            collection = destination$iv$iv2;
            boolean bl = false;
            aminoAcidFragment = it.toAminoAcidFragment();
            collection.add(aminoAcidFragment);
        }
        List unfilteredAminoAcidFragments = (List) destination$iv$iv2;
        Iterable $receiver$iv3 = unfilteredAminoAcidFragments;
        destination$iv$iv2 = $receiver$iv3;
        Collection destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it3;
            AminoAcidFragment bl = (AminoAcidFragment) item$iv$iv;
            collection = destination$iv$iv3;
            boolean bl2 = false;
            aminoAcidFragment = this.enrich((AminoAcidFragment) it3, highQualityAminoAcidCounts);
            collection.add(aminoAcidFragment);
        }
        List result = (List) destination$iv$iv3;
        return result;

         */
    }

    private final AminoAcidFragment enrich(AminoAcidFragment fragment, SequenceCount count)
    {
        return fragment;

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        List<Integer> initialIntersect = fragment.aminoAcidLoci();
        Iterable iterable = $receiver$iv = (Iterable) initialIntersect;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            String actual;
            int loci = ((Number) element$iv$iv).intValue();
            boolean bl = false;
            Collection<String> allowed = count.sequenceAt(loci);
            if(!allowed.contains(actual = fragment.aminoAcid(loci)))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List filteredIntersect = (List) destination$iv$iv;
        return fragment.intersectAminoAcidLoci(filteredIntersect);
         */
    }

    public AminoAcidQualEnrichment(int minEvidence)
    {
        this.mMinEvidence = minEvidence;
    }
}
