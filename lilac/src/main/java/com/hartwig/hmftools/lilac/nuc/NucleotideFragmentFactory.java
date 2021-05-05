package com.hartwig.hmftools.lilac.nuc;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.read.AminoAcidIndices;
import com.hartwig.hmftools.lilac.sam.SAMCodingRecord;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.LociPosition;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class NucleotideFragmentFactory
{
    private final Map<HlaSequenceLoci,SuffixTree> mInsertSuffixTrees;
    private final Map<HlaSequenceLoci,SuffixTree> mDeleteSuffixTrees;
    private final int mMinBaseQuality;
    private final LociPosition mLociPosition;

    public NucleotideFragmentFactory(
            int minBaseQuality, final List<HlaSequenceLoci> inserts, final List<HlaSequenceLoci> deletes, final LociPosition lociPosition)
    {
        mMinBaseQuality = minBaseQuality;
        mLociPosition = lociPosition;

        mInsertSuffixTrees = Maps.newHashMap();
        mDeleteSuffixTrees = Maps.newHashMap();

        inserts.stream().forEach(x -> mInsertSuffixTrees.put(x, new SuffixTree(x.sequence())));
        deletes.stream().forEach(x -> mDeleteSuffixTrees.put(x, new SuffixTree(x.sequence())));
    }

    public final NucleotideFragment createAlignmentFragments(final SAMCodingRecord samCoding, @NotNull NamedBed codingRegion)
    {
        return null;

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;

        Iterable iterable = $receiver$iv = (Iterable) samCoding.alignmentsOnly();
        Collection destination$iv$iv = new ArrayList();
        NucleotideFragment $receiver$iv$iv$iv = $receiver$iv$iv;
        Iterator iterator = $receiver$iv$iv$iv.iterator();
        while(iterator.hasNext())
        {
            NucleotideFragment nucleotideFragment;
            Object element$iv$iv$iv;
            Object element$iv$iv = element$iv$iv$iv = iterator.next();
            com.hartwig.hmftools.lilac.sam.SAMCodingRecord it = (com.hartwig.hmftools.lilac.sam.SAMCodingRecord) element$iv$iv;
            boolean bl = false;
            if(createFragment(it, codingRegion) == null)
            {
                continue;
            }
            NucleotideFragment it$iv$iv = nucleotideFragment;
            destination$iv$iv.add(it$iv$iv);
        }
        List all = (List) destination$iv$iv;
        if(all.isEmpty())
        {
            return null;
        }
        $receiver$iv = all;
        Iterator iterator$iv = $receiver$iv.iterator();
        if(!iterator$iv.hasNext())
        {
            throw (Throwable) new UnsupportedOperationException("Empty collection can't be reduced.");
        }
        Object accumulator$iv = iterator$iv.next();
        while(iterator$iv.hasNext())
        {
            void y;
            $receiver$iv$iv$iv = (NucleotideFragment) iterator$iv.next();
            NucleotideFragment x = (NucleotideFragment) accumulator$iv;
            boolean bl = false;
            accumulator$iv = NucleotideFragment.Companion.merge(x, (NucleotideFragment) y);
        }
        return (NucleotideFragment) accumulator$iv;

         */
    }

    public final NucleotideFragment createFragment(@NotNull SAMCodingRecord samCoding, @NotNull NamedBed codingRegion)
    {
        return null;

        /*
        void $receiver$iv$iv;
        char[] $receiver$iv;
        String string;
        Collection collection;
        Object aminoAcids;
        int n;
        int nucleotideStartLoci;
        Intrinsics.checkParameterIsNotNull((Object) samCoding, (String) "samCoding");
        Intrinsics.checkParameterIsNotNull((Object) codingRegion, (String) "codingRegion");
        boolean reverseStrand = samCoding.getReverseStrand();
        int samCodingStartLoci = reverseStrand
                ? mLociPosition.nucelotideLoci(samCoding.getPositionEnd())
                : mLociPosition.nucelotideLoci(samCoding.getPositionStart());
        int samCodingEndLoci = reverseStrand
                ? mLociPosition.nucelotideLoci(samCoding.getPositionStart())
                : mLociPosition.nucelotideLoci(samCoding.getPositionEnd());
        char[] codingRegionRead = samCoding.codingRegionRead(reverseStrand);
        int[] codingRegionQuality = samCoding.codingRegionQuality(reverseStrand);
        if(samCoding.containsIndel() || samCoding.containsSoftClip())
        {
            IntRange aminoAcidIndices = AminoAcidIndices.INSTANCE.indices(samCodingStartLoci, samCodingEndLoci);
            nucleotideStartLoci = aminoAcidIndices.getFirst() * 3;
            String sequence2 =
                    ArraysKt.joinToString$default((char[]) codingRegionRead, (CharSequence) "", null, null, (int) 0, null, null, (int) 62, null);
            CharSequence charSequence = sequence2;
            n = nucleotideStartLoci - samCodingStartLoci;
            String string2 = charSequence;
            if(string2 == null)
            {
                throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
            }
            String string3 = string2.substring(n);
            Intrinsics.checkExpressionValueIsNotNull((Object) string3, (String) "(this as java.lang.String).substring(startIndex)");
            aminoAcids = Codons.aminoAcids((String) string3);
            Intrinsics.checkExpressionValueIsNotNull((Object) aminoAcids, (String) "aminoAcids");
            charSequence = (CharSequence) aminoAcids;
            if(charSequence.length() > 0)
            {
                Iterable $receiver$iv$iv2;
                Iterable $receiver$iv2;
                Object result;
                Object element$iv$iv3;
                Collection destination$iv$iv;
                Pair pair;
                Object it;
                Object item$iv$iv;
                Iterable $receiver$iv$iv3;
                Iterable $receiver$iv3;
                n = aminoAcidIndices.getFirst() - samCoding.getSoftClippedStart() / 3 - samCoding.maxIndelSize();
                IntRange matchRangeAllowed =
                        new IntRange(n, aminoAcidIndices.getFirst() + samCoding.maxIndelSize() + samCoding.getSoftClippedEnd() / 3);
                Map<HlaSequenceLoci, SuffixTree> map = mInsertSuffixTrees;
                void var15_18 = $receiver$iv3;
                Collection destination$iv$iv2 = new ArrayList($receiver$iv3.size());
                Iterator iterator = $receiver$iv$iv3;
                Iterator iterator2 = iterator.entrySet().iterator();
                while(iterator2.hasNext())
                {
                    void it2;
                    Pair pair2 = item$iv$iv = iterator2.next();
                    collection = destination$iv$iv2;
                    boolean bl2 = false;
                    string = new Pair(it2.getKey(), (Object) ((SuffixTree) it2.getValue()).indices((String) aminoAcids));
                    collection.add(string);
                }
                $receiver$iv3 = (List) destination$iv$iv2;
                $receiver$iv$iv3 = $receiver$iv3;
                destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
                for(Object item$iv$iv2 : $receiver$iv$iv3)
                {
                    void $receiver$iv$iv4;
                    void $receiver$iv4;
                    Pair pair3;
                    item$iv$iv = (Pair) item$iv$iv2;
                    collection = destination$iv$iv2;
                    boolean bl = false;
                    Object object = it.getFirst();
                    Object object2 = it.getSecond();
                    Intrinsics.checkExpressionValueIsNotNull((Object) object2, (String) "it.second");
                    Iterable bl2 = (Iterable) object2;
                    Object $i$f$mapTo = object;
                    Pair $i$f$map = pair3;
                    Pair pair4 = pair3;
                    pair = $receiver$iv4;
                    destination$iv$iv = new ArrayList();
                    for(Object element$iv$iv2 : $receiver$iv$iv4)
                    {
                        Integer i = (Integer) element$iv$iv2;
                        boolean bl3 = false;
                        Integer n2 = i;
                        Intrinsics.checkExpressionValueIsNotNull((Object) n2, (String) "i");
                        if(!matchRangeAllowed.contains(n2.intValue()))
                        {
                            continue;
                        }
                        destination$iv$iv.add(element$iv$iv2);
                    }
                    List list = (List) destination$iv$iv;
                    $i$f$map($i$f$mapTo, (Object) list);
                    string = pair4;
                    collection.add(string);
                }
                $receiver$iv3 = (List) destination$iv$iv2;
                $receiver$iv$iv3 = $receiver$iv3;
                destination$iv$iv2 = new ArrayList();
                for(Object element$iv$iv3 : $receiver$iv$iv3)
                {
                    it = (Pair) element$iv$iv3;
                    boolean bl = false;
                    Collection $receiver$iv4 = (Collection) it.getSecond();
                    if(!(!$receiver$iv4.isEmpty()))
                    {
                        continue;
                    }
                    destination$iv$iv2.add(element$iv$iv3);
                }
                List matchingInserts = (List) destination$iv$iv2;
                $receiver$iv3 = matchingInserts;
                if(!$receiver$iv3.isEmpty())
                {
                    Pair best = (Pair) matchingInserts.get(0);
                    String string4 = samCoding.getId();
                    Object e = ((List) best.getSecond()).get(0);
                    Intrinsics.checkExpressionValueIsNotNull(e, (String) "best.second[0]");
                    result =
                            createNucleotideSequence(string4, codingRegion, ((Number) e).intValue(), (String) aminoAcids, (HlaSequenceLoci) best
                                    .getFirst());
                    return result;
                }
                result = mDeleteSuffixTrees;
                destination$iv$iv2 = $receiver$iv2;
                Collection destination$iv$iv3 = new ArrayList($receiver$iv2.size());
                element$iv$iv3 = $receiver$iv$iv2;
                it = element$iv$iv3.entrySet().iterator();
                while(it.hasNext())
                {
                    void it3;
                    Map.Entry item$iv$iv3;
                    Map.Entry $receiver$iv4 = item$iv$iv3 = (Map.Entry) it.next();
                    collection = destination$iv$iv3;
                    boolean bl4 = false;
                    string = new Pair(it3.getKey(), (Object) ((SuffixTree) it3.getValue()).indices((String) aminoAcids));
                    collection.add(string);
                }
                $receiver$iv2 = (List) destination$iv$iv3;
                $receiver$iv$iv2 = $receiver$iv2;
                destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
                element$iv$iv3 = $receiver$iv$iv2.iterator();
                while(element$iv$iv3.hasNext())
                {
                    void $receiver$iv$iv5;
                    void $receiver$iv5;
                    Pair pair5;
                    void it4;
                    Object item$iv$iv3 = item$iv$iv = element$iv$iv3.next();
                    collection = destination$iv$iv3;
                    boolean bl = false;
                    Object object = it4.getFirst();
                    Object object3 = it4.getSecond();
                    Intrinsics.checkExpressionValueIsNotNull((Object) object3, (String) "it.second");
                    Iterable bl4 = (Iterable) object3;
                    Object $i$f$mapTo = object;
                    Pair $i$f$map = pair5;
                    pair = pair5;
                    destination$iv$iv = $receiver$iv5;
                    Collection destination$iv$iv4 = new ArrayList();
                    for(Object element$iv$iv4 : $receiver$iv$iv5)
                    {
                        Integer i = (Integer) element$iv$iv4;
                        boolean bl5 = false;
                        Integer n3 = i;
                        Intrinsics.checkExpressionValueIsNotNull((Object) n3, (String) "i");
                        if(!matchRangeAllowed.contains(n3.intValue()))
                        {
                            continue;
                        }
                        destination$iv$iv4.add(element$iv$iv4);
                    }
                    List list = (List) destination$iv$iv4;
                    $i$f$map($i$f$mapTo, (Object) list);
                    string = pair;
                    collection.add(string);
                }
                $receiver$iv2 = (List) destination$iv$iv3;
                $receiver$iv$iv2 = $receiver$iv2;
                destination$iv$iv3 = new ArrayList();
                for(Object element$iv$iv5 : $receiver$iv$iv2)
                {
                    Pair it5 = (Pair) element$iv$iv5;
                    boolean bl = false;
                    Collection collection2 = (Collection) it5.getSecond();
                    if(!(!collection2.isEmpty()))
                    {
                        continue;
                    }
                    destination$iv$iv3.add(element$iv$iv5);
                }
                List matchingDeletes = (List) destination$iv$iv3;
                $receiver$iv2 = matchingDeletes;
                if(!$receiver$iv2.isEmpty())
                {
                    Pair best = (Pair) matchingDeletes.get(0);
                    String string5 = samCoding.getId();
                    Object e = ((List) best.getSecond()).get(0);
                    Intrinsics.checkExpressionValueIsNotNull(e, (String) "best.second[0]");
                    NucleotideFragment result2 =
                            createNucleotideSequence(string5, codingRegion, ((Number) e).intValue(), (String) aminoAcids, (HlaSequenceLoci) best
                                    .getFirst());
                    return result2;
                }
            }
            if(samCoding.containsIndel())
            {
                return null;
            }
        }
        if(samCodingStartLoci < 0 || samCodingEndLoci < 0)
        {
            return null;
        }
        nucleotideStartLoci = samCodingStartLoci;
        List loci = CollectionsKt.toList((Iterable) ((Iterable) new IntRange(nucleotideStartLoci, samCodingEndLoci)));
        aminoAcids = $receiver$iv = codingRegionRead;
        Collection destination$iv$iv = new ArrayList($receiver$iv.length);
        n = ((void) $receiver$iv$iv).length;
        for(int i = 0; i < n; ++i)
        {
            void it;
            void item$iv$iv;
            void result2 = item$iv$iv = $receiver$iv$iv[i];
            collection = destination$iv$iv;
            boolean bl = false;
            string = String.valueOf((char) it);
            collection.add(string);
        }
        List nucleotides2 = (List) destination$iv$iv;
        List qualities = ArraysKt.toList((int[]) codingRegionQuality);
        return new NucleotideFragment(samCoding.getId(), SetsKt.setOf((Object) codingRegion.name()), loci, qualities, nucleotides2);

         */
    }

    private final NucleotideFragment createNucleotideSequence(
            String id, NamedBed codingRegion, int startLoci, String bamSequence, HlaSequenceLoci hlaSequence)
    {
        return null;

        /*
        void $receiver$iv$iv;
        Object object;
        Collection collection;
        Iterable $receiver$iv$iv2;
        int n;
        void $receiver$iv$iv3;
        Iterable $receiver$iv;
        int endLoci = endLoci(startLoci, bamSequence, hlaSequence);
        int n2 = startLoci;
        List aminoAcidLoci = CollectionsKt.toList((Iterable) ((Iterable) new IntRange(n2, endLoci)));
        Iterable iterable = $receiver$iv = (Iterable) aminoAcidLoci;
        Iterable destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv3)
        {
            int it = ((Number) element$iv$iv).intValue();
            n = 0;
            Iterable list$iv$iv = CollectionsKt.listOf((Object[]) new Integer[] { 3 * it, 3 * it + 1, 3 * it + 2 });
            CollectionsKt.addAll((Collection) destination$iv$iv, (Iterable) list$iv$iv);
        }
        List nucleotideLoci = (List) destination$iv$iv;
        Iterable $receiver$iv2 = aminoAcidLoci;
        destination$iv$iv = $receiver$iv2;
        Iterable destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv2)
        {
            void it;
            n = ((Number) item$iv$iv).intValue();
            collection = destination$iv$iv2;
            boolean bl = false;
            object = hlaSequence.sequence((int) it);
            collection.add(object);
        }
        $receiver$iv2 = (List) destination$iv$iv2;
        $receiver$iv$iv2 = $receiver$iv2;
        destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            String it = (String) element$iv$iv;
            boolean bl = false;
            Iterable list$iv$iv = Companion.createNucleotidesFromAminoAcid(it);
            CollectionsKt.addAll((Collection) destination$iv$iv2, (Iterable) list$iv$iv);
        }
        List nucleotides2 = (List) destination$iv$iv2;
        Iterable $receiver$iv3 = nucleotideLoci;
        destination$iv$iv2 = $receiver$iv3;
        Collection destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            int bl = ((Number) item$iv$iv).intValue();
            collection = destination$iv$iv3;
            boolean bl2 = false;
            object = mMinBaseQuality;
            collection.add(object);
        }
        List qualities = (List) destination$iv$iv3;
        Set genes = SetsKt.setOf((Object) codingRegion.name());
        return new NucleotideFragment(id, genes, nucleotideLoci, qualities, nucleotides2);

         */
    }

    private final int endLoci(int startLoci, String bamSequence, HlaSequenceLoci hlaSequence)
    {
        return 0;

        /*
        StringBuilder builder = new StringBuilder();
        int n = startLoci;
        int n2 = hlaSequence.length();
        while(n < n2)
        {
            void loci;
            builder.append(hlaSequence.getSequences().get((int) loci));
            String string = builder.toString();
            Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "builder.toString()");
            if(!StringsKt.startsWith$default((String) bamSequence, (String) string, (boolean) false, (int) 2, null))
            {
                return (int) (loci - true);
            }
            ++loci;
        }
        return hlaSequence.length() - 1;

         */
    }

    public static List<String> createNucleotidesFromAminoAcid(@NotNull String aminoAcid)
    {
        return Lists.newArrayList();

        /*
        Intrinsics.checkParameterIsNotNull((Object) aminoAcid, (String) "aminoAcid");
        if(Intrinsics.areEqual((Object) aminoAcid, (Object) "."))
        {
            return CollectionsKt.listOf((Object[]) new String[] { ".", ".", "." });
        }
        String codons = Codons.codons((String) aminoAcid);
        Object[] objectArray = new String[3];
        objectArray[0] = String.valueOf(codons.charAt(0));
        objectArray[1] = String.valueOf(codons.charAt(1));
        String string = codons;
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "codons");
        String string2 = string;
        int n = 2;
        String string3 = string2;
        if(string3 == null)
        {
            throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
        }
        String string4 = string3.substring(n);
        Intrinsics.checkExpressionValueIsNotNull((Object) string4, (String) "(this as java.lang.String).substring(startIndex)");
        objectArray[2] = string4;
        return CollectionsKt.listOf((Object[]) objectArray);

         */
    }
}
