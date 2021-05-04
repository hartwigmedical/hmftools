package com.hartwig.hmftools.lilac.evidence;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public final class PhasedEvidence implements Comparable<PhasedEvidence>
{
    private final int[] mAminoAcidIndices;
    private final Map<String,Integer> mEvidenceMap;

    public PhasedEvidence(final int[] aminoAcidIndices, final Map<String, Integer> evidence)
    {
        mAminoAcidIndices = aminoAcidIndices;
        mEvidenceMap = evidence;
    }

    public final int[] getAminoAcidIndices() { return mAminoAcidIndices; }
    public final Map<String, Integer> getEvidence()
    {
        return mEvidenceMap;
    }

    private boolean consistentWithAny(final List<HlaSequenceLoci> candidates, final String sequence)
    {
        return candidates.stream().anyMatch(x -> x.consistentWith(sequence, mAminoAcidIndices));
    }

    public PhasedEvidence inconsistentEvidence(final List<HlaSequenceLoci> candidates)
    {
        /*
            fun consistentWithAny(sequence: String): Boolean {
            return candidates.any { it.consistentWith(sequence, *aminoAcidIndices) }
        }
        return PhasedEvidence(aminoAcidIndices, evidence.filter { !consistentWithAny(it.key) })

         */

        // TODO
        PhasedEvidence phasedEvidence = null;

        return phasedEvidence;
    }

    public final int[] unambiguousHeadIndices()
    {
        return null;
        /*
        void var3_3;
        void $receiver$iv$iv;
        int[] $receiver$iv;
        int[] nArray = $receiver$iv = mAminoAcidIndices;
        Collection destination$iv$iv = new ArrayList();
        void $receiver$iv$iv$iv = $receiver$iv$iv;
        int index$iv$iv$iv = 0;
        for(void item$iv$iv$iv : $receiver$iv$iv$iv)
        {
            void element$iv$iv;
            int n = index$iv$iv$iv++;
            void var9_9 = item$iv$iv$iv;
            int index$iv$iv = n;
            void var11_11 = element$iv$iv;
            int index = index$iv$iv;
            boolean bl = false;
            if(!(index < unambiguousHeadLength()))
            {
                continue;
            }
            destination$iv$iv.add((int) element$iv$iv);
        }
        return CollectionsKt.toIntArray((Collection) ((List) var3_3));

         */
    }

    public final int[] unambiguousTailIndices()
    {
        return null;

        /*
        void $receiver$iv$iv;
        int[] $receiver$iv;
        int minIndex = mAminoAcidIndices.length - unambiguousTailLength();
        int[] nArray = $receiver$iv = mAminoAcidIndices;
        Collection destination$iv$iv = new ArrayList();
        void $receiver$iv$iv$iv = $receiver$iv$iv;
        int index$iv$iv$iv = 0;
        for(void item$iv$iv$iv : $receiver$iv$iv$iv)
        {
            void element$iv$iv;
            int n = index$iv$iv$iv++;
            void var10_10 = item$iv$iv$iv;
            int index$iv$iv = n;
            void var12_12 = element$iv$iv;
            int index = index$iv$iv;
            boolean bl = false;
            if(!(index >= minIndex))
            {
                continue;
            }
            destination$iv$iv.add((int) element$iv$iv);
        }
        return CollectionsKt.toIntArray((Collection) ((List) destination$iv$iv));

         */
    }

    private final int unambiguousTailLength()
    {
        int endIndex = mAminoAcidIndices.length;

        for(int i = 0; i < mAminoAcidIndices.length; ++i)
        {
            int startIndex = endIndex - i - 1;

            Set<String> evidenceTails = mEvidenceMap.keySet().stream().map(x -> x.substring(startIndex, endIndex)).collect(Collectors.toSet());

            if(evidenceTails.size() == mEvidenceMap.size())
                return i + 1;
        }

        return mAminoAcidIndices.length;
    }

    private final int unambiguousHeadLength()
    {
        for(int i = 0; i < mAminoAcidIndices.length; ++i)
        {
            int length = i + 1;
            Set<String> evidenceTails = mEvidenceMap.keySet().stream().map(x -> x.substring(0, length)).collect(Collectors.toSet());

            if(evidenceTails.size() == mEvidenceMap.size())
                return length;
        }

        return mAminoAcidIndices.length;
    }

    public final boolean contains(final PhasedEvidence other)
    {
        // TODO
        /*
        Intrinsics.checkParameterIsNotNull((Object) other, (String) "other");
        int overlap =
                CollectionsKt.intersect((Iterable) ArraysKt.toSet((int[]) other.mAminoAcidIndices), (Iterable) ArraysKt.toSet((int[]) mAminoAcidIndices))
                        .size();
        return overlap == other.mAminoAcidIndices.length;

         */

        return false;
    }

    public final int minEvidence()
    {
        return mEvidenceMap.values().stream().mapToInt(x -> x).min().orElse(0);
    }

    public final int totalEvidence()
    {
        return mEvidenceMap.values().stream().mapToInt(x -> x).sum();
    }

    public final PhasedEvidence removeSingles(int minTotalEvidence)
    {
        if(totalEvidence() >= minTotalEvidence)
        {
            Map<String, Integer> newEvidence = Maps.newHashMap();
            mEvidenceMap.entrySet().stream().filter(x -> x.getValue() > 1).forEach(x -> newEvidence.put(x.getKey(), x.getValue()));
            return new PhasedEvidence(mAminoAcidIndices, newEvidence);
        }

        return this;
    }

    public final String evidenceString()
    {
        StringJoiner sj = new StringJoiner(" ");

        mEvidenceMap.entrySet().forEach(x -> sj.add(toEvidenceString(x.getKey(), x.getValue())));

        return sj.toString();
    }

    public final String toEvidenceString(final String evidence, int count)
    {
        /*
                val resultBuilder = StringJoiner("")
        var evidenceIndex = 0
        for (i in aminoAcidIndices[0]..aminoAcidIndices[aminoAcidIndices.lastIndex]) {
            if (i in aminoAcidIndices) {
                resultBuilder.add(evidence[evidenceIndex].toString())
                evidenceIndex++
            } else {
                resultBuilder.add("-")
            }
        }

        return resultBuilder.add("=").add(count.toString()).toString()

         */
        // TODO
        return "";
    }

    public String toString()
    {
        int uniqueTail = unambiguousTailLength();
        int uniqueHead = unambiguousHeadLength();
        StringBuilder stringBuilder = new StringBuilder().append("PhasedEvidence(head=")
                .append(uniqueHead)
                .append(" tail=")
                .append(uniqueTail)
                .append(" loci=")
                .append(mAminoAcidIndices.length)
                .append(" types=")
                .append(mEvidenceMap.size())
                .append(" indices=");
        int[] nArray = mAminoAcidIndices;
        String string = Arrays.toString(nArray);

        return stringBuilder.append(string)
                .append(", evidence=")
                .append(mEvidenceMap)
                .append(", total=")
                .append(totalEvidence())
                .append(')')
                .toString();
    }

    @Override
    public int compareTo(final PhasedEvidence other)
    {
        // TODO - careful
        /*
        int totalCompare = -Intrinsics.compare((int) totalEvidence(), (int) other.totalEvidence());
        if(totalCompare != 0)
        {
            return totalCompare;
        }
        return Intrinsics.compare((int) minEvidence(), (int) other.minEvidence());

         */

        return 1;
    }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof PhasedEvidence))
            return false;

        // TODO check
        int[] nArray = mAminoAcidIndices;
        int[] nArray2 = ((PhasedEvidence) other).mAminoAcidIndices;
        return Arrays.equals(nArray, nArray2);
    }

    public int hashCode()
    {
        int[] nArray = mAminoAcidIndices;
        return Arrays.hashCode(nArray);
    }

    public static PhasedEvidence evidence(final List<AminoAcidFragment> aminoAcidFragments, final int[] indices)
    {
        // TODO:
        return null;

        /*
                    val filteredFragments = aminoAcidFragments.filter { it.containsAll(indices) }
            val aminoAcidEvidence = filteredFragments.map { it.toAminoAcids(indices) }.groupingBy { it }.eachCount()
            return PhasedEvidence(indices, aminoAcidEvidence)

         */
    }

    private static String toAminoAcids(final AminoAcidFragment fragment, final int[] indices)
    {
        // TODO
        return "";
        // return indices.map { this.aminoAcid(it) }.joinToString("")

    }

    private static boolean containsAll(AminoAcidFragment $receiver, int[] indices)
    {
        // TODO:
        return false;
        // return indices.all { this.containsAminoAcid(it) }
    }
}
