package com.hartwig.hmftools.lilac.qc;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilackt.SequenceCount;
import com.hartwig.hmftools.lilackt.qc.Haplotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class HaplotypeQC
{
    private final int mUnusedHaplotypes;
    private final int mUnusedHaplotypeMaxSupport;
    private final int mUnusedHaplotypeMaxLength;
    private final int mUnusedHaplotypesPon;

    private static final int MIN_SUPPORT = 7;

    private static final List<Haplotype> PON_HAPLOTYPES = Lists.newArrayList();

    @NotNull
    public final List<String> header()
    {
        return Lists.newArrayList(
                "unusedHaplotypes", "unusedHaplotypeMaxSupport", "unusedHaplotypeMaxLength", "unusedHaplotypesPon");
    }

    @NotNull
    public final List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(mUnusedHaplotypes),
                String.valueOf(mUnusedHaplotypeMaxSupport), 
                String.valueOf(mUnusedHaplotypeMaxLength),
                String.valueOf(mUnusedHaplotypesPon));
    }

    public final int getUnusedHaplotypes()
    {
        return mUnusedHaplotypes;
    }

    public HaplotypeQC(int unusedHaplotypes, int unusedHaplotypeMaxSupport, int unusedHaplotypeMaxLength, int unusedHaplotypesPon)
    {
        mUnusedHaplotypes = unusedHaplotypes;
        mUnusedHaplotypeMaxSupport = unusedHaplotypeMaxSupport;
        mUnusedHaplotypeMaxLength = unusedHaplotypeMaxLength;
        mUnusedHaplotypesPon = unusedHaplotypesPon;

        /* TODO - populate PON_HAPLOTYPES
        HaplotypeQC.class.getResource("/pon/haplotypes.txt")
            .readText()
            .split("\n")
            .map { Haplotype.fromString(it) }
         */

    }

    private final boolean inPon(Haplotype $receiver)
    {
        boolean bl;
        block3:
        {
            Iterable $receiver$iv = PON_HAPLOTYPES;
            if($receiver$iv instanceof Collection && ((Collection) $receiver$iv).isEmpty())
            {
                bl = false;
            }
            else
            {
                for(Object element$iv : $receiver$iv)
                {
                    Haplotype it = (Haplotype) element$iv;
                    boolean bl2 = false;
                    if(!it.contains($receiver))
                    {
                        continue;
                    }
                    bl = true;
                    break block3;
                }
                bl = false;
            }
        }
        return bl;
    }

    public static HaplotypeQC create(
            int minEvidence, final Set<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> winners, final List<com.hartwig.hmftools.lilac.evidence.PhasedEvidence> evidence,
            final SequenceCount aminoAcidCount)
    {
        return null;

        /*
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;

        Iterable iterable = evidence;
        void var7_7 = $receiver$iv;
        Object destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            com.hartwig.hmftools.lilac.evidence.PhasedEvidence it = (com.hartwig.hmftools.lilac.evidence.PhasedEvidence) element$iv$iv;
            boolean bl = false;
            Iterable list$iv$iv = Companion.unmatchedHaplotype(it, minEvidence, (Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci>) winners, aminoAcidCount);
            CollectionsKt.addAll((Collection) destination$iv$iv, (Iterable) list$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                Haplotype it = (Haplotype) b;
                boolean bl = false;
                Comparable comparable = Integer.valueOf(it.getSupportingFragments());
                it = (Haplotype) a;
                Comparable comparable2 = comparable;
                bl = false;
                Integer n = it.getSupportingFragments();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) n);
            }
        };
        $receiver$iv = CollectionsKt.sortedWith((Iterable) $receiver$iv$iv, (Comparator) destination$iv$iv);
        Iterable<String> set$iv = new HashSet();
        Object list$iv = new ArrayList();
        for(Object e$iv : $receiver$iv)
        {
            Haplotype x = (Haplotype) e$iv;
            boolean bl = false;
            String key$iv = String.valueOf(x.getStartLocus()) + String.valueOf(x.getEndLocus()) + x.getHaplotype();
            if(!((HashSet) set$iv).add(key$iv))
            {
                continue;
            }
            ((ArrayList) list$iv).add(e$iv);
        }
        $receiver$iv = (List) list$iv;
        set$iv = $receiver$iv;
        list$iv = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                Haplotype it = (Haplotype) a;
                boolean bl = false;
                Comparable comparable = Integer.valueOf(it.getStartLocus());
                it = (Haplotype) b;
                Comparable comparable2 = comparable;
                bl = false;
                Integer n = it.getStartLocus();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) n);
            }
        };
        List allUnmatched = CollectionsKt.sortedWith(set$iv, (Comparator) list$iv);
        int pon = 0;
        int unusedCount = 0;
        int maxSupport = 0;
        int maxLength = 0;
        for(Haplotype unmatched : allUnmatched)
        {
            if(unmatched.getSupportingFragments() < 7)
            {
                continue;
            }
            if(inPon(unmatched))
            {
                ++pon;
                logger.info("    UNMATCHED_PON_HAPLTOYPE - " + unmatched);
                continue;
            }
            int n = unmatched.getSupportingFragments();
            maxSupport = Math.max(maxSupport, n);
            n = unmatched.getHaplotype().length();
            maxLength = Math.max(maxLength, n);
            ++unusedCount;
            logger.warn("    UNMATCHED_HAPLTOYPE - " + unmatched);
        }
        return new HaplotypeQC(unusedCount, maxSupport, maxLength, pon);

         */
    }

    public final List<Haplotype> unmatchedHaplotype(@NotNull com.hartwig.hmftools.lilac.evidence.PhasedEvidence $receiver, int minEvidence,
            @NotNull Collection<com.hartwig.hmftools.lilac.seq.HlaSequenceLoci> winners, @NotNull SequenceCount aminoAcidCount)
    {
        return Lists.newArrayList();

        /*
        Object object;
        Object object2;
        boolean bl;
        Map.Entry it;
        Object $receiver$iv$iv;
        Object $receiver$iv;
        Intrinsics.checkParameterIsNotNull((Object) $receiver, (String) "$receiver");
        Intrinsics.checkParameterIsNotNull(winners, (String) "winners");
        Intrinsics.checkParameterIsNotNull((Object) aminoAcidCount, (String) "aminoAcidCount");
        Function1<String, Boolean> $fun$consistentWithAny$1 = new Function1<String, Boolean>($receiver, winners)
        {
            public final boolean invoke(@NotNull String sequence2)
            {
                boolean bl;
                block3:
                {
                    Intrinsics.checkParameterIsNotNull((Object) sequence2, (String) "sequence");
                    Iterable $receiver$iv = $winners;
                    if($receiver$iv instanceof Collection && ((Collection) $receiver$iv).isEmpty())
                    {
                        bl = false;
                    }
                    else
                    {
                        for(T element$iv : $receiver$iv)
                        {
                            com.hartwig.hmftools.lilac.seq.HlaSequenceLoci it = (HlaSequenceLoci) element$iv;
                            boolean bl2 = false;
                            int[] nArray = receiver$0.getAminoAcidIndices();
                            if(!it.consistentWith(sequence2, Arrays.copyOf(nArray, nArray.length)))
                            {
                                continue;
                            }
                            bl = true;
                            break block3;
                        }
                        bl = false;
                    }
                }
                return bl;
            }

            {
                receiver$0 = phasedEvidence;
                $winners = collection;
                super(1);
            }
        };
        Map<String, Integer> map = $receiver.getEvidence();
        void var8_7 = $receiver$iv;
        Object destination$iv$iv = new LinkedHashMap();
        Object object3 = $receiver$iv$iv;
        Iterator iterator = object3.entrySet().iterator();
        while(iterator.hasNext())
        {
            Map.Entry entry;
            it = entry = iterator.next();
            bl = false;
            if(!(!$fun$consistentWithAny$1.invoke((String) it.getKey())))
            {
                continue;
            }
            destination$iv$iv.put(entry.getKey(), entry.getValue());
        }
        $receiver$iv = destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new LinkedHashMap();
        object3 = $receiver$iv$iv;
        for(Map.Entry entry : object3.entrySet())
        {
            it = entry;
            bl = false;
            if(!(((Number) it.getValue()).intValue() >= minEvidence))
            {
                continue;
            }
            destination$iv$iv.put(entry.getKey(), entry.getValue());
        }
        Object unmatched = destination$iv$iv;
        if(unmatched.isEmpty())
        {
            return CollectionsKt.emptyList();
        }
        $receiver$iv$iv = $receiver$iv = unmatched;
        destination$iv$iv = new ArrayList($receiver$iv.size());
        object3 = $receiver$iv$iv;
        for(Map.Entry entry : object3.entrySet())
        {
            it = entry;
            object2 = destination$iv$iv;
            bl = false;
            object = new Pair(it.getKey(), it.getValue());
            object2.add(object);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        object3 = $receiver$iv$iv.iterator();
        while(object3.hasNext())
        {
            Object item$iv$iv2 = object3.next();
            Pair pair = (Pair) item$iv$iv2;
            object2 = destination$iv$iv;
            boolean bl2 = false;
            object = Haplotype.Companion.create($receiver.getAminoAcidIndices(), (Pair<String, Integer>) pair, aminoAcidCount);
            object2.add(object);
        }
        return (List) destination$iv$iv;
        */
    }
}
