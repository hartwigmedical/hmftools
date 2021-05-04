package com.hartwig.hmftools.lilac;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class SequenceCount
{
    // private final int length;
    private final int mMinCount;
    private final Map<String,Integer>[] mSeqCountsList;

    public SequenceCount(int minCount, final Map<String, Integer>[] seqCounts)
    {
        mMinCount = minCount;
        mSeqCountsList = seqCounts;
        // length = mSeqCounts.length;
    }

    public final int getLength()
    {
        return mSeqCountsList.length;
    }

    public final Map<String,Integer> get(int locus) { return mSeqCountsList[locus]; }

    public final List<Integer> heterozygousLoci()
    {
        List<Integer> indices = Lists.newArrayList();

        for(int i = 0; i < mSeqCountsList.length; ++i)
        {
            if(isHeterozygous(i))
                indices.add(i);
        }

        return indices;
    }

    public final List<Integer> homozygousIndices()
    {
        List<Integer> indices = Lists.newArrayList();

        for(int i = 0; i < mSeqCountsList.length; ++i)
        {
            if(isHomozygous(i))
                indices.add(i);
        }

        return indices;
    }

    private final boolean isHomozygous(int index)
    {
        Map<String,Integer> seqCounts = mSeqCountsList[index];
        return seqCounts.values().stream().filter(x -> x >= mMinCount).count() == 1;
    }

    private final boolean isHeterozygous(int index)
    {
        Map<String,Integer> seqCounts = mSeqCountsList[index];
        return seqCounts.values().stream().filter(x -> x >= mMinCount).count() > 1;
    }

    public final List<String> sequenceAt(int index)
    {
        if (index >= mSeqCountsList.length)
            return Lists.newArrayList();

        return mSeqCountsList[index].entrySet().stream()
                .filter(x -> x.getValue() >= mMinCount)
                .map(x -> x.getKey()).collect(Collectors.toList());
    }

    public final int depth(int index)
    {
        Map<String,Integer> seqCounts = mSeqCountsList[index];
        return seqCounts.values().stream().mapToInt(x -> x).sum();
    }

    public final void writeVertically(final String fileName)
    {
        // TODO
        /*
        Intrinsics.checkParameterIsNotNull((Object) fileName, (String) "fileName");
        File file = new File(fileName);
        FilesKt.writeText$default((File) file, (String) "", null, (int) 2, null);
        int n = 0;
        int n2 = length;
        while(n < n2)
        {
            void $receiver$iv$iv;
            Object $receiver$iv;
            void i;
            StringJoiner lineBuilder = new StringJoiner("\t").add(String.valueOf((int) i));
            Object object = $receiver$iv = mSeqCountsList[i];
            Object destination$iv$iv322 = new ArrayList($receiver$iv.size());
            void var10_15 = $receiver$iv$iv;
            Pair pair = var10_15.entrySet().iterator();
            while(pair.hasNext())
            {
                void k;
                void $dstr$k$v;
                Map.Entry item$iv$iv;
                Map.Entry entry = item$iv$iv = pair.next();
                Collection collection = destination$iv$iv322;
                boolean bl = false;
                void var15_20 = $dstr$k$v;
                String string = (String) var15_20.getKey();
                var15_20 = $dstr$k$v;
                int v = ((Number) var15_20.getValue()).intValue();
                Pair pair2 = new Pair((Object) k, (Object) v);
                collection.add(pair2);
            }
            $receiver$iv = (List) destination$iv$iv322;
            int $receiver$iv2 = 0;
            int destination$iv$iv322 = 5;
            object = $receiver$iv;
            destination$iv$iv322 = new Comparator<T>()
            {

                public final int compare(T a, T b)
                {
                    Pair it = (Pair) a;
                    boolean bl = false;
                    Comparable comparable = Integer.valueOf(((Number) it.getSecond()).intValue());
                    it = (Pair) b;
                    Comparable comparable2 = comparable;
                    bl = false;
                    Integer n = ((Number) it.getSecond()).intValue();
                    return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) n);
                }
            };
            List baseCountList =
                    CollectionsKt.reversed((Iterable) CollectionsKt.sortedWith((Iterable) object, (Comparator) destination$iv$iv322));
            int $i$f$sortedBy = baseCountList.size() - 1;
            int n3 = Math.min(destination$iv$iv322, $i$f$sortedBy);
            if($receiver$iv2 <= n3)
            {
                void j;
                do
                {
                    void base;
                    pair = (Pair) baseCountList.get((int) (++j));
                    String destination$iv$iv322 = (String) pair.component1();
                    int count = ((Number) pair.component2()).intValue();
                    lineBuilder.add((CharSequence) base).add(String.valueOf(count));
                } while(j != n3);
            }
            FilesKt.appendText$default((File) file, (String) (lineBuilder.toString() + "\n"), null, (int) 2, null);
            ++i;
        }
         */
    }

    private static void increment(Map<String,Integer>[] seqCountsList, int index, String aminoAcid)
    {
        Map<String,Integer> seqCounts = seqCountsList[index];

        // CHECK
        Integer count = seqCounts.get(aminoAcid);
        if(count != null)
            seqCounts.put(aminoAcid, count + 1);
        else
            seqCounts.put(aminoAcid, 1);

        // seqCounts.compute(aminoAcid, 1)
    }

    public static SequenceCount nucleotides(int minCount, final List<NucleotideFragment> fragments)
    {
        int length = fragments.stream().mapToInt(x -> x.maxLoci()).max().orElse(0) + 1;

        Map<String,Integer>[] seqCountsList = new Map[length];

        for(com.hartwig.hmftools.lilac.nuc.NucleotideFragment fragment : fragments)
        {
            for(Integer lociIndex : fragment.getNucleotideLoci())
            {
                String nucleotide = fragment.getNucleotides().get(lociIndex);
                increment(seqCountsList, lociIndex, nucleotide);
            }
        }

        return new SequenceCount(minCount, seqCountsList);
    }

    public static SequenceCount aminoAcids(int minCount, final List<AminoAcidFragment> fragments)
    {
        int length = fragments.stream().mapToInt(x -> x.maxLoci()).max().orElse(0) + 1;

        Map<String,Integer>[] seqCountsList = new Map[length];

        for(NucleotideFragment fragment : fragments)
        {
            for(Integer lociIndex : fragment.getNucleotideLoci())
            {
                String nucleotide = fragment.getNucleotides().get(lociIndex);
                increment(seqCountsList, lociIndex, nucleotide);
            }
        }

        return new SequenceCount(minCount, seqCountsList);

        /*
            val length = fragments.map { it.aminoAcidLoci().max() ?: -1 }.max()!! + 1
            val count = Array(length) { mutableMapOf<String, Int>() }

            for (fragment in fragments) {
                for (index in fragment.aminoAcidLoci()) {
                    val aminoAcid = fragment.aminoAcid(index)
                    count.increment(index, aminoAcid)
                }
            }
            return SequenceCount(minCount, Array(length) { count[it] })

         */

    }

}
