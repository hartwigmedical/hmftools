package com.hartwig.hmftools.lilac.evidence;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public final class PhasedEvidence implements Comparable<PhasedEvidence>
{
    private final int[] mAminoAcidIndices;
    private final List<Integer> mAminoAcidIndexList;
    private final Map<String,Integer> mEvidenceMap;

    public PhasedEvidence(final int[] aminoAcidIndices, final Map<String,Integer> evidence)
    {
        mAminoAcidIndices = aminoAcidIndices;
        mAminoAcidIndexList = Lists.newArrayList();
        Arrays.stream(aminoAcidIndices).forEach(x -> mAminoAcidIndexList.add(x));
        mEvidenceMap = evidence;
    }

    public PhasedEvidence(final List<Integer> aminoAcidIndexList, final Map<String,Integer> evidence)
    {
        mAminoAcidIndexList = aminoAcidIndexList;
        mAminoAcidIndices = new int[aminoAcidIndexList.size()];

        for(int i = 0; i < aminoAcidIndexList.size(); ++i)
            mAminoAcidIndices[i] = aminoAcidIndexList.get(i);

        mEvidenceMap = evidence;
    }

    public final int[] getAminoAcidIndices() { return mAminoAcidIndices; }
    public final List<Integer> getAminoAcidIndexList() { return mAminoAcidIndexList; }

    public final Map<String,Integer> getEvidence()
    {
        return mEvidenceMap;
    }

    private boolean consistentWithAny(final List<HlaSequenceLoci> candidates, final String sequence)
    {
        return candidates.stream().anyMatch(x -> x.consistentWith(sequence, mAminoAcidIndices));
    }

    public PhasedEvidence inconsistentEvidence(final List<HlaSequenceLoci> candidates)
    {
        Map<String,Integer> filteredEvidence = mEvidenceMap.entrySet().stream()
                .filter(x -> !consistentWithAny(candidates, x.getKey()))
                .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

        return new PhasedEvidence(mAminoAcidIndices, filteredEvidence);
    }

    public List<Integer> unambiguousHeadIndices()
    {
        int headLength = unambiguousHeadLength();

        List<Integer> indices = Lists.newArrayList();

        for(int i = 0; i < min(mAminoAcidIndices.length, headLength); ++i)
        {
            indices.add(mAminoAcidIndices[i]);
        }

        return indices;
    }

    public List<Integer> unambiguousTailIndices()
    {
        int minIndex = mAminoAcidIndices.length - unambiguousTailLength();

        List<Integer> indices = Lists.newArrayList();

        for(int i = minIndex; i < mAminoAcidIndices.length; ++i)
        {
            indices.add(mAminoAcidIndices[i]);
        }

        return indices;
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
        Set<Integer> overlapIndices = Sets.newHashSet();
        Arrays.stream(mAminoAcidIndices).forEach(x -> overlapIndices.add(x));
        Arrays.stream(other.getAminoAcidIndices()).forEach(x -> overlapIndices.add(x));

        int overlapCount = overlapIndices.size();
        return overlapCount == other.getAminoAcidIndices().length;
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
            Map<String, Integer> newEvidence = mEvidenceMap.entrySet().stream()
                    .filter(x -> x.getValue() > 1).collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

            return new PhasedEvidence(mAminoAcidIndices, newEvidence);
        }

        return this;
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
        int totalCompare = other.totalEvidence() - totalEvidence();
        if(totalCompare != 0)
            return totalCompare;

        int minCompare = minEvidence() - other.minEvidence();

        if(minCompare != 0)
            return minCompare > 0 ? 1 : -1;

        return 0;
    }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof PhasedEvidence))
            return false;

        final int[] otherIndices = ((PhasedEvidence)other).getAminoAcidIndices();

        if(mAminoAcidIndices.length != otherIndices.length)
            return false;

        for(int i = 0; i < mAminoAcidIndices.length; ++i)
        {
            if(mAminoAcidIndices[i] != otherIndices[i])
                return false;
        }

        return true;
    }

    public int hashCode()
    {
        return Arrays.hashCode(mAminoAcidIndices);
    }

    public static PhasedEvidence evidence(final List<Fragment> fragments, final List<Integer> indices)
    {
        List<Fragment> filteredFragments = fragments.stream().filter(x -> x.containsAminoAcids(indices)).collect(Collectors.toList());
        final Map<String,Integer> evidence = Maps.newHashMap();

        for(Fragment fragment : filteredFragments)
        {
            String sequence = fragment.aminoAcids(indices);
            Integer count = evidence.get(sequence);
            evidence.put(sequence, count != null ? count + 1 : 1);
        }

        return new PhasedEvidence(indices, evidence);
    }

    public static void logInconsistentEvidence(final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        List<HlaSequenceLoci> expectedSequences = candidates.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

        for (HlaSequenceLoci sequence : expectedSequences)
        {
            for (PhasedEvidence phasedEvidence : evidence)
            {
                if (!sequence.consistentWith(phasedEvidence))
                {
                    LL_LOGGER.warn("Expected allele {} filtered by {}", sequence.Allele, phasedEvidence);
                }
            }
        }
    }

}
