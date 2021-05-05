package com.hartwig.hmftools.lilac.seq;

import static com.hartwig.hmftools.lilac.LilacConstants.DEL_STR;
import static com.hartwig.hmftools.lilac.LilacConstants.WILD_CHAR;
import static com.hartwig.hmftools.lilac.LilacConstants.WILD_STR;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class HlaSequenceLoci
{
    private final HlaAllele mAllele;
    private final List<String> mSequences;

    public HlaSequenceLoci(final HlaAllele allele, final List<String> sequences)
    {
        mAllele = allele;
        mSequences = Lists.newArrayList();
        mSequences.addAll(sequences);
    }

    public int length() { return mSequences.size(); }

    public final HlaAllele getAllele() { return mAllele;}

    public final List<String> getSequences() { return mSequences; }

    public final boolean containsInserts() { return mSequences.stream().anyMatch(x -> x.length() > 1); }

    public final boolean containsDeletes() { return mSequences.stream().anyMatch(x -> x.equals(DEL_STR)); }

    public final boolean containsIndels() { return containsDeletes() || containsInserts(); }

    public final String sequence(int locus) { return mSequences.get(locus); }

    public final String sequence(int startLocus, int endLocus)
    {
        // CHECK
        StringJoiner sj = new StringJoiner("");
        for(int i = startLocus; i <= endLocus; ++i)
        {
            sj.add(mSequences.get(i));
        }

        return sj.toString();
    }

    public final String sequence(final Set<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");

        for(Integer index : indices)
        {
            if(index < mSequences.size())
                sj.add(mSequences.get(index));
            else
                sj.add(WILD_STR);
        }

        return sj.toString();
    }

    public final String sequence()
    {
        StringJoiner sj = new StringJoiner("");
        for(String sequence : mSequences)
        {
            sj.add(sequence.replace(DEL_STR, ""));
        }

        return sj.toString();
    }

    public String toString()
    {
        return String.format("HlaSequenceLoci(allele=%s, sequence=%s)", mAllele, sequence());
    }

    public final boolean consistentWith(final PhasedEvidence evidence)
    {
        return consistentWithAny(
                evidence.getEvidence().keySet().stream().collect(Collectors.toList()), evidence.getAminoAcidIndices());
    }

    public final boolean consistentWithAny(final Set<String> targetSequences, final int[] targetIndices)
    {
        return consistentWithAny(targetSequences.stream().collect(Collectors.toList()), targetIndices);
    }

    public final boolean consistentWithAny(final List<String> targetSequences, final int[] targetIndices)
    {
        Set<Integer> tiSet = Sets.newHashSet();
        Arrays.stream(targetIndices).forEach(x -> tiSet.add(x));
        return consistentWithAny(targetSequences, tiSet);
    }

    public final boolean consistentWithAny(final List<String> targetSequences, final Set<Integer> targetIndices)
    {
        return targetSequences.stream().anyMatch(x -> match(x, targetIndices) != HlaSequenceMatch.NONE);
    }

    public final boolean consistentWith(final String targetSequence, final Set<Integer> targetIndices)
    {
        return match(targetSequence, targetIndices) != HlaSequenceMatch.NONE;
    }

    // TODO - require both methods?
    public final boolean consistentWith(final String targetSequence, final int[] targetIndices)
    {
        Set<Integer> tiSet = Sets.newHashSet();
        Arrays.stream(targetIndices).forEach(x -> tiSet.add(x));
        return match(targetSequence, tiSet) != HlaSequenceMatch.NONE;
    }

    public final HlaSequenceMatch match(final String targetSequence, final Set<Integer> targetIndices)
    {
        if(targetIndices.isEmpty())
            return HlaSequenceMatch.NONE;

        String hlaSequence = sequence(targetIndices);

        if(hlaSequence.length() != targetSequence.length())
            return HlaSequenceMatch.NONE;

        int wildCardCount = 0;

        for(int i = 0; i < targetSequence.length(); ++i)
        {
            char targetChr = targetSequence.charAt(i);
            char hlaSeqChr = hlaSequence.charAt(i);

            if(hlaSeqChr != WILD_CHAR && hlaSeqChr != targetChr)
                return HlaSequenceMatch.NONE;

            if (hlaSeqChr == WILD_CHAR)
                wildCardCount++;
        }

        if (wildCardCount > 0)
            return wildCardCount == targetIndices.size() ? HlaSequenceMatch.WILD : HlaSequenceMatch.PARTIAL;

        return HlaSequenceMatch.FULL;
    }

    public static List<HlaSequenceLoci> create(final List<HlaSequence> sequences)
    {
        final String reference = sequences.get(0).getRawSequence();

        return sequences.stream()
                .map(x -> create(x.Allele, x.getRawSequence(), reference))
                .filter(x -> !x.getSequences().isEmpty())
                .collect(Collectors.toList());
    }

    public static HlaSequenceLoci create(final HlaAllele allele, final String sequence, final String reference)
    {
        final List<String> sequences = Lists.newArrayList();

        int insLength = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            char seqChar = sequence.charAt(i);
            boolean isBaseInserted = seqChar != '.' && (i >= reference.length() || reference.charAt(i) == '.');
            boolean isBaseIgnored = (seqChar == '.' && i < reference.length() && reference.charAt(i) == '.') || (seqChar == '|');

            if(insLength > 0 && !isBaseInserted)
            {
                String insertStr = sequence.substring(i - insLength, i);
                sequences.set(sequences.size() - 1, sequences.get(sequences.size() - 1) + insertStr);
                //sequences[sequences.size - 1] = sequences[sequences.size - 1] + insertStr
                insLength = 0;
            }

            if(isBaseInserted)
            {
                insLength++;
            }
            else if(!isBaseIgnored)
            {
                // String locusSequence = seqChar; // sequence[i].toString()

                if(seqChar == '-')
                    sequences.add(String.valueOf(reference.charAt(i)));
                else
                    sequences.add(String.valueOf(seqChar));
            }
        }

        return new HlaSequenceLoci(allele, sequences);
    }
}
