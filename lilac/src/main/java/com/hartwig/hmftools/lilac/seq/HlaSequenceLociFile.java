package com.hartwig.hmftools.lilac.seq;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci;

public class HlaSequenceLociFile
{
    // writes exon boundaries, and then the allele sequences
    /*
        HLA-C Boundaries                                |                                                                                         |                                                                                           |
        A*01:01                 MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAA
     */

    public void write(
            final String file, final Set<Integer> aBoundaries, final Set<Integer> bBoundaries,
            final Set<Integer> cBoundaries, final List<HlaSequenceLoci> sequences)
    {
        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        boolean bl;

        Collection collection = sequences;
        boolean bl2 = bl = !collection.isEmpty();
        if(!bl)
        {
            String string = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(string.toString());
        }
        File outputFile = new File(file);
        FilesKt.writeText$default((File) outputFile, (String) "", null, (int) 2, null);
        Iterable iterable = $receiver$iv = (Iterable) sequences;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            HlaSequenceLoci hlaSequenceLoci = (HlaSequenceLoci) item$iv$iv;
            Collection collection2 = destination$iv$iv;
            boolean bl3 = false;
            List<Integer> list = INSTANCE.lengths((HlaSequenceLoci) it);
            collection2.add(list);
        }
        $receiver$iv = (List) destination$iv$iv;
        Iterator iterator$iv22 = $receiver$iv.iterator();
        if(!iterator$iv22.hasNext())
        {
            throw (Throwable) new UnsupportedOperationException("Empty collection can't be reduced.");
        }
        Object accumulator$iv = iterator$iv22.next();
        while(iterator$iv22.hasNext())
        {
            void right;
            List list = (List) iterator$iv22.next();
            List left = (List) accumulator$iv;
            boolean bl4 = false;
            accumulator$iv = INSTANCE.maxLengths(left, (List<Integer>) right);
        }
        List maxLengths = (List) accumulator$iv;
        String templateSequence = this.padInserts(sequences.get(0), maxLengths);
        FilesKt.appendText$default((File) outputFile, (String) (StringsKt.padEnd((String) "HLA-A Boundaries", (int) 20, (char) ' ') + "\t"
                + this.boundaryString((Collection<Integer>) aBoundaries, maxLengths) + "\n"), null, (int) 2, null);
        FilesKt.appendText$default((File) outputFile, (String) (StringsKt.padEnd((String) "HLA-B Boundaries", (int) 20, (char) ' ') + "\t"
                + this.boundaryString((Collection<Integer>) bBoundaries, maxLengths) + "\n"), null, (int) 2, null);
        FilesKt.appendText$default((File) outputFile, (String) (StringsKt.padEnd((String) "HLA-C Boundaries", (int) 20, (char) ' ') + "\t"
                + this.boundaryString((Collection<Integer>) cBoundaries, maxLengths) + "\n"), null, (int) 2, null);
        FilesKt.appendText$default((File) outputFile, (String) (
                StringsKt.padEnd((String) String.valueOf(sequences.get(0).getAllele()), (int) 20, (char) ' ') + "\t" + templateSequence
                        + "\n"), null, (int) 2, null);
        int iterator$iv22 = 1;
        int n = sequences.size();
        while(iterator$iv22 < n)
        {
            void i;
            String victimSequence = this.diff(this.padInserts(sequences.get((int) i), maxLengths), templateSequence);
            FilesKt.appendText$default((File) outputFile, (String) (
                    StringsKt.padEnd((String) String.valueOf(sequences.get((int) i).getAllele()), (int) 20, (char) ' ') + "\t"
                            + victimSequence + "\n"), null, (int) 2, null);
            ++i;
        }

         */
    }

    private final String boundaryString(Collection<Integer> boundaries, List<Integer> lengths)
    {
        /*
            val range = (0..boundaries.max()!!)
            val unpadded = range.mapIndexed { index, _ -> if (index in boundaries) "|" else " " }
            return unpadded.mapIndexed { index, value -> value.padEnd(lengths[index], ' ') }.joinToString("")
         */

        return "";

        /*
        String string;
        int index;
        Collection collection;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        boolean bl;
        Collection<Integer> collection2 = boundaries;
        boolean bl2 = bl = !collection2.isEmpty();
        if(!bl)
        {
            String string2 = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(string2.toString());
        }
        int n = 0;
        Comparable comparable = CollectionsKt.max((Iterable) boundaries);
        if(comparable == null)
        {
            Intrinsics.throwNpe();
        }
        IntRange range = new IntRange(n, ((Number) ((Object) comparable)).intValue());
        Iterable iterable = $receiver$iv = (Iterable) range;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        int index$iv$iv = 0;
        Iterator iterator = $receiver$iv$iv.iterator();
        while(iterator.hasNext())
        {
            int item$iv$iv = ((IntIterator) iterator).nextInt();
            int n2 = index$iv$iv++;
            int n3 = item$iv$iv;
            int n4 = n2;
            collection = destination$iv$iv;
            boolean bl3 = false;
            string = boundaries.contains(index) ? "|" : " ";
            collection.add(string);
        }
        List unpadded = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv = (Iterable) unpadded;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        index$iv$iv = 0;
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void value;
            int n5 = index$iv$iv++;
            String $noName_1 = (String) item$iv$iv;
            index = n5;
            collection = destination$iv$iv;
            boolean bl4 = false;
            string = StringsKt.padEnd((String) value, (int) ((Number) lengths.get(index)).intValue(), (char) ' ');
            collection.add(string);
        }
        return CollectionsKt.joinToString$default((Iterable) ((List) destination$iv$iv), (CharSequence) "", null, null, (int) 0, null, null, (int) 62, null);
        */
    }

    /*
    private final String diff(String victim, String reference)
    {
        void $receiver$iv$iv;
        CharSequence $receiver$iv = victim;
        CharSequence charSequence = $receiver$iv;
        Collection destination$iv$iv = new ArrayList($receiver$iv.length());
        int index$iv$iv = 0;
        void var7_7 = $receiver$iv$iv;
        for(int i = 0; i < var7_7.length(); ++i)
        {
            void index;
            char item$iv$iv = var7_7.charAt(i);
            int n = index$iv$iv++;
            char c = item$iv$iv;
            int n2 = n;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            Character c2 = Character.valueOf(INSTANCE.diff(victim.charAt((int) index),
                    index < reference.length() ? reference.charAt((int) index) : (char) '!'));
            collection.add(c2);
        }
        return CollectionsKt.joinToString$default((Iterable) ((List) destination$iv$iv), (CharSequence) "", null, null, (int) 0, null, null, (int) 62, null);
    }

    private final char diff(char victim, char reference)
    {
        char c = victim;
        return (char) (c == '|' ? 124 : (c == '*' ? 42 : (c == '.' ? 46 : (c == reference ? 45 : (int) victim))));
    }

    private final String padInserts(HlaSequenceLoci $receiver, List<Integer> lengths)
    {
        CharSequence charSequence;
        block3:
        {
            CharSequence $receiver$iv$iv;
            Object $receiver$iv = $receiver.getSequences();
            Iterable iterable = $receiver$iv;
            Collection destination$iv$iv22 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
            int index$iv$iv = 0;
            Iterator iterator = $receiver$iv$iv.iterator();
            while(iterator.hasNext())
            {
                void index;
                void value;
                Object item$iv$iv = iterator.next();
                int n = index$iv$iv++;
                String string = (String) item$iv$iv;
                int n2 = n;
                Collection collection = destination$iv$iv22;
                boolean bl = false;
                String string2 = StringsKt.padEnd((String) value, (int) ((Number) lengths.get((int) index)).intValue(), (char) '.');
                collection.add(string2);
            }
            $receiver$iv =
                    CollectionsKt.joinToString$default((Iterable) ((List) destination$iv$iv22), (CharSequence) "", null, null, (int) 0, null, null, (int) 62, null);
            $receiver$iv$iv = (CharSequence) $receiver$iv;
            int destination$iv$iv22 = $receiver$iv$iv.length();
            boolean bl = false;
            while(--destination$iv$iv22 >= 0)
            {
                void index$iv$iv2;
                char x = $receiver$iv$iv.charAt((int) index$iv$iv2);
                boolean bl2 = false;
                if(!(x == '.'))
                {
                    charSequence = $receiver$iv$iv.subSequence(0, (int) (index$iv$iv2 + true));
                    break block3;
                }
                --index$iv$iv2;
            }
            charSequence = "";
        }
        return ((Object) charSequence).toString();
    }

    private final List<Integer> lengths(HlaSequenceLoci $receiver)
    {
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) $receiver.getSequences();
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            String string = (String) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            Integer n = it.length();
            collection.add(n);
        }
        return (List) destination$iv$iv;
    }

    private final List<Integer> maxLengths(List<Integer> left, List<Integer> right)
    {
        List result = new ArrayList();
        int n = 0;
        int n2 = left.size();
        int n3 = right.size();
        int n4 = Math.max(n2, n3);
        while(n < n4)
        {
            void i;
            int leftValue = i < left.size() ? ((Number) left.get((int) i)).intValue() : 0;
            int rightValue = i < right.size() ? ((Number) right.get((int) i)).intValue() : 0;
            result.add(Math.max(leftValue, rightValue));
            ++i;
        }
        return result;
    }
    */
}
