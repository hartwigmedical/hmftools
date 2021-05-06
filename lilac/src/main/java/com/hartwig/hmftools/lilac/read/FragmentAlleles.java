package com.hartwig.hmftools.lilac.read;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.Collection;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class FragmentAlleles
{
    private final AminoAcidFragment mFragment;
    private final List<HlaAllele> mFull;
    private final List<HlaAllele> mPartial;
    private final List<HlaAllele> mWild;
    
    public final boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mPartial.contains(allele) || mWild.contains(allele);
    }

    public final AminoAcidFragment getFragment() { return mFragment; }

    public final Collection<HlaAllele> getFull() { return mFull; }

    public final Collection<HlaAllele> getPartial() { return mPartial; }

    public final Collection<HlaAllele> getWild() { return mWild; }

    public FragmentAlleles(
            final AminoAcidFragment fragment, final List<HlaAllele> full, final List<HlaAllele> partial, final List<HlaAllele> wild)
    {
        mFragment = fragment;
        mFull = full;
        mPartial = partial;
        mWild = wild;
    }

    private static FragmentAlleles filter(final FragmentAlleles $receiver, final List<HlaAllele> alleles)
    {
        return null;

        /*
        HlaAllele it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        FragmentAlleles fragmentAlleles;
        Iterable iterable = $receiver.getFull();
        AminoAcidFragment aminoAcidFragment = $receiver.getFragment();
        FragmentAlleles fragmentAlleles2 = fragmentAlleles;
        FragmentAlleles fragmentAlleles3 = fragmentAlleles;
        void var4_7 = $receiver$iv;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!alleles.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection = (List) destination$iv$iv;
        $receiver$iv = $receiver.getPartial();
        collection = collection;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!alleles.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection2 = (List) destination$iv$iv;
        $receiver$iv = $receiver.getWild();
        collection2 = collection2;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!alleles.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List list = (List) destination$iv$iv;
        fragmentAlleles2(aminoAcidFragment, collection, collection2, list);
        return fragmentAlleles3;

         */
    }

    public final List<FragmentAlleles> filter(@NotNull List<FragmentAlleles> $receiver, @NotNull Collection<HlaAllele> alleles)
    {
        return Lists.newArrayList();

        /*
        FragmentAlleles it;
        Iterable $receiver$iv$iv;
        Intrinsics.checkParameterIsNotNull($receiver, (String) "$receiver");
        Intrinsics.checkParameterIsNotNull(alleles, (String) "alleles");
        Iterable $receiver$iv = $receiver;
        Iterable iterable = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            FragmentAlleles fragmentAlleles = (FragmentAlleles) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            FragmentAlleles fragmentAlleles2 = Companion.filter(it, alleles);
            collection.add(fragmentAlleles2);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (FragmentAlleles) element$iv$iv;
            boolean bl = false;
            Collection<HlaAllele> collection = it.getFull();
            if(!(!collection.isEmpty() || !(collection = it.getPartial()).isEmpty()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        return (List) destination$iv$iv;

         */
    }

    public static List<FragmentAlleles> create(final List<AminoAcidFragment> aminoAcidFragments, final Collection<Integer> hetLoci,
            final Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> sequences, final Collection<Integer> nucleotideLoci,
            final Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> nucleotideSequences)
    {
        return Lists.newArrayList();

        /*
        FragmentAlleles it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull(aminoAcidFragments, (String) "aminoAcidFragments");
        Intrinsics.checkParameterIsNotNull(hetLoci, (String) "hetLoci");
        Intrinsics.checkParameterIsNotNull(sequences, (String) "sequences");
        Intrinsics.checkParameterIsNotNull(nucleotideLoci, (String) "nucleotideLoci");
        Intrinsics.checkParameterIsNotNull(nucleotideSequences, (String) "nucleotideSequences");
        Iterable iterable = aminoAcidFragments;
        void var7_7 = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            AminoAcidFragment aminoAcidFragment = (AminoAcidFragment) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            FragmentAlleles fragmentAlleles =
                    Companion.create((AminoAcidFragment) ((Object) it), hetLoci, sequences, nucleotideLoci, nucleotideSequences);
            collection.add(fragmentAlleles);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (FragmentAlleles) element$iv$iv;
            boolean bl = false;
            Collection<HlaAllele> collection = it.getFull();
            if(!(!collection.isEmpty() || !(collection = it.getPartial()).isEmpty()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        return (List) destination$iv$iv;

         */
    }

    private static FragmentAlleles create(AminoAcidFragment aminoAcidFragment, Collection<Integer> aminoAcidLoci,
            Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> aminoAcidSequences, Collection<Integer> nucleotideLoci,
            Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> nucleotideSequences)
    {
        return null;

        /*
        void $receiver$iv$iv;
        void $receiver$iv$iv2;
        void $receiver$iv$iv3;
        Iterable $receiver$iv$iv4;
        Iterable $receiver$iv$iv5;
        Iterable $receiver$iv$iv6;
        Pair it;
        Iterable $receiver$iv$iv7;
        Iterable $receiver$iv;
        Iterable $receiver$iv$iv8;
        Iterable $receiver$iv2;
        Object item$iv$iv3;
        Iterable $receiver$iv$iv9;
        Object object;
        Pair it2;
        Collection collection;
        Iterable $receiver$iv$iv10;
        Iterable $receiver$iv3;
        int[] fragmentNucleotideLoci =
                CollectionsKt.toIntArray((Collection) CollectionsKt.sorted((Iterable) CollectionsKt.intersect((Iterable) aminoAcidFragment
                        .getNucleotideLoci(), (Iterable) nucleotideLoci)));
        String fragmentNucleotides =
                aminoAcidFragment.nucleotides(Arrays.copyOf(fragmentNucleotideLoci, fragmentNucleotideLoci.length));
        Iterable iterable = nucleotideSequences;
        void var10_9 = $receiver$iv3;
        Iterable destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
        for(Object item$iv$iv2 : $receiver$iv$iv10)
        {
            com.hartwig.hmftools.lilac.seq.HlaSequenceLoci hlaSequenceLoci = (com.hartwig.hmftools.lilac.seq.HlaSequenceLoci) item$iv$iv2;
            collection = destination$iv$iv;
            boolean bl = false;
            object =
                    new Pair((Object) it2.getAllele(), (Object) it2.match(fragmentNucleotides, Arrays.copyOf(fragmentNucleotideLoci, fragmentNucleotideLoci.length)));
            collection.add(object);
        }
        $receiver$iv3 = (List) destination$iv$iv;
        $receiver$iv$iv10 = $receiver$iv3;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv10)
        {
            it2 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it2.getSecond()) != HlaSequenceMatch.NONE))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv3 = (List) destination$iv$iv;
        $receiver$iv$iv10 = $receiver$iv3;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv10)
        {
            it2 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!aminoAcidFragment.getGenes().contains("HLA-" + ((HlaAllele) it2.getFirst()).getGene()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv3 = (List) destination$iv$iv;
        $receiver$iv$iv10 = $receiver$iv3;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
        for(Object item$iv$iv2 : $receiver$iv$iv10)
        {
            it2 = (Pair) item$iv$iv2;
            collection = destination$iv$iv;
            boolean bl = false;
            object = new Pair((Object) ((HlaAllele) it2.getFirst()).asFourDigit(), it2.getSecond());
            collection.add(object);
        }
        List matchingNucleotideSequences = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv));
        Iterable $receiver$iv4 = matchingNucleotideSequences;
        destination$iv$iv = $receiver$iv4;
        Collection destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv9)
        {
            Pair it3 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it3.getSecond()) == HlaSequenceMatch.FULL))
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        $receiver$iv4 = (List) destination$iv$iv2;
        $receiver$iv$iv9 = $receiver$iv4;
        destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv4, (int) 10));
        for(Object item$iv$iv3 : $receiver$iv$iv9)
        {
            Pair it3 = (Pair) item$iv$iv3;
            collection = destination$iv$iv2;
            boolean bl = false;
            object = (HlaAllele) it3.getFirst();
            collection.add(object);
        }
        Set fullNucleotideMatch = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv2));
        $receiver$iv$iv9 = matchingNucleotideSequences;
        destination$iv$iv2 = $receiver$iv2;
        Collection destination$iv$iv3 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv8)
        {
            Pair it4 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it4.getSecond()) == HlaSequenceMatch.PARTIAL
                    || (HlaSequenceMatch) ((Object) it4.getSecond()) == HlaSequenceMatch.WILD))
            {
                continue;
            }
            destination$iv$iv3.add(element$iv$iv);
        }
        $receiver$iv2 = (List) destination$iv$iv3;
        $receiver$iv$iv8 = $receiver$iv2;
        destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
        for(Object item$iv$iv4 : $receiver$iv$iv8)
        {
            Pair it4 = (Pair) item$iv$iv4;
            collection = destination$iv$iv3;
            boolean bl = false;
            object = (HlaAllele) it4.getFirst();
            collection.add(object);
        }
        Set partialNucleotideMatch =
                CollectionsKt.subtract((Iterable) CollectionsKt.toSet((Iterable) ((List) destination$iv$iv3)), (Iterable) fullNucleotideMatch);
        int[] fragmentAminoAcidLoci =
                CollectionsKt.toIntArray((Collection) CollectionsKt.sorted((Iterable) CollectionsKt.intersect((Iterable) aminoAcidFragment
                        .aminoAcidLoci(), (Iterable) aminoAcidLoci)));
        String fragmentAminoAcids = aminoAcidFragment.aminoAcids(Arrays.copyOf(fragmentAminoAcidLoci, fragmentAminoAcidLoci.length));
        item$iv$iv3 = aminoAcidSequences;
        void item$iv$iv4 = $receiver$iv;
        Iterable<Object> destination$iv$iv4 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv5 : $receiver$iv$iv7)
        {
            com.hartwig.hmftools.lilac.seq.HlaSequenceLoci $i$f$map = (HlaSequenceLoci) item$iv$iv5;
            collection = destination$iv$iv4;
            boolean bl = false;
            object =
                    new Pair((Object) it.getAllele(), (Object) it.match(fragmentAminoAcids, Arrays.copyOf(fragmentAminoAcidLoci, fragmentAminoAcidLoci.length)));
            collection.add(object);
        }
        $receiver$iv = (List) destination$iv$iv4;
        $receiver$iv$iv7 = $receiver$iv;
        destination$iv$iv4 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv7)
        {
            it = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it.getSecond()) != HlaSequenceMatch.NONE))
            {
                continue;
            }
            destination$iv$iv4.add(element$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv4;
        $receiver$iv$iv7 = $receiver$iv;
        destination$iv$iv4 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv7)
        {
            it = (Pair) element$iv$iv;
            boolean bl = false;
            if(!aminoAcidFragment.getGenes().contains("HLA-" + ((HlaAllele) it.getFirst()).getGene()))
            {
                continue;
            }
            destination$iv$iv4.add(element$iv$iv);
        }
        List matchingAminoAcidSequences = (List) destination$iv$iv4;
        Iterable $receiver$iv5 = matchingAminoAcidSequences;
        destination$iv$iv4 = $receiver$iv5;
        Iterable destination$iv$iv5 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv6)
        {
            Pair it5 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it5.getSecond()) == HlaSequenceMatch.FULL))
            {
                continue;
            }
            destination$iv$iv5.add(element$iv$iv);
        }
        $receiver$iv5 = (List) destination$iv$iv5;
        $receiver$iv$iv6 = $receiver$iv5;
        destination$iv$iv5 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv5, (int) 10));
        for(Object item$iv$iv6 : $receiver$iv$iv6)
        {
            Pair it5 = (Pair) item$iv$iv6;
            collection = destination$iv$iv5;
            boolean bl = false;
            object = (HlaAllele) it5.getFirst();
            collection.add(object);
        }
        Set fullAminoAcidMatch = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv5));
        Iterable $receiver$iv6 = matchingAminoAcidSequences;
        destination$iv$iv5 = $receiver$iv6;
        Iterable destination$iv$iv6 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv5)
        {
            Pair it6 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it6.getSecond()) == HlaSequenceMatch.PARTIAL))
            {
                continue;
            }
            destination$iv$iv6.add(element$iv$iv);
        }
        $receiver$iv6 = (List) destination$iv$iv6;
        $receiver$iv$iv5 = $receiver$iv6;
        destination$iv$iv6 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv6, (int) 10));
        for(Object item$iv$iv7 : $receiver$iv$iv5)
        {
            Pair it6 = (Pair) item$iv$iv7;
            collection = destination$iv$iv6;
            boolean bl = false;
            object = (HlaAllele) it6.getFirst();
            collection.add(object);
        }
        Set partialAminoAcidMatch = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv6));
        Iterable $receiver$iv7 = matchingAminoAcidSequences;
        destination$iv$iv6 = $receiver$iv7;
        Iterable destination$iv$iv7 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv4)
        {
            Pair it7 = (Pair) element$iv$iv;
            boolean bl = false;
            if(!((HlaSequenceMatch) ((Object) it7.getSecond()) == HlaSequenceMatch.WILD))
            {
                continue;
            }
            destination$iv$iv7.add(element$iv$iv);
        }
        $receiver$iv7 = (List) destination$iv$iv7;
        $receiver$iv$iv4 = $receiver$iv7;
        destination$iv$iv7 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv7, (int) 10));
        for(Object item$iv$iv8 : $receiver$iv$iv4)
        {
            Pair it7 = (Pair) item$iv$iv8;
            collection = destination$iv$iv7;
            boolean bl = false;
            object = (HlaAllele) it7.getFirst();
            collection.add(object);
        }
        Set wildAminoAcidMatch = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv7));
        if(fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty())
        {
            return new FragmentAlleles(aminoAcidFragment, fullAminoAcidMatch, partialAminoAcidMatch, wildAminoAcidMatch);
        }
        Iterable $receiver$iv8 = fullAminoAcidMatch;
        destination$iv$iv7 = $receiver$iv8;
        Iterable destination$iv$iv8 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv3)
        {
            HlaAllele it8 = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!fullNucleotideMatch.contains(it8.asFourDigit()))
            {
                continue;
            }
            destination$iv$iv8.add(element$iv$iv);
        }
        List consistentFull = (List) destination$iv$iv8;
        Iterable $receiver$iv9 = fullAminoAcidMatch;
        destination$iv$iv8 = $receiver$iv9;
        Iterable destination$iv$iv9 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            HlaAllele it9 = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!partialNucleotideMatch.contains(it9.asFourDigit()))
            {
                continue;
            }
            destination$iv$iv9.add(element$iv$iv);
        }
        List downgradedToPartial = (List) destination$iv$iv9;
        Iterable $receiver$iv10 = partialAminoAcidMatch;
        destination$iv$iv9 = $receiver$iv10;
        Collection destination$iv$iv10 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            HlaAllele it10 = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!partialNucleotideMatch.contains(it10.asFourDigit()))
            {
                continue;
            }
            destination$iv$iv10.add(element$iv$iv);
        }
        List otherPartial = (List) destination$iv$iv10;
        return new FragmentAlleles(aminoAcidFragment, consistentFull, CollectionsKt.union((Iterable) downgradedToPartial, (Iterable) otherPartial), wildAminoAcidMatch);

         */
    }
}
