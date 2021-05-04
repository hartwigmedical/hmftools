package com.hartwig.hmftools.lilac.seq;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.Arrays;
import java.util.Collection;
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

    public final boolean containsInserts() { return mSequences.stream().anyMatch(x -> x.length() > 0); }

    public final boolean containsDeletes() { return mSequences.stream().anyMatch(x -> x.equals(".")); }

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
        for(int i = 0; i <= mSequences.size(); ++i)
        {
            if(indices.contains(i))
                sj.add(mSequences.get(i));
        }

        return sj.toString();
    }

    public final String sequence()
    {
        // TODO
        StringJoiner sj = new StringJoiner("");
        for(int i = 0; i <= mSequences.size(); ++i)
        {
            sj.add(mSequences.get(i).replace(".", ""));
        }

        return sj.toString();
    }

    public String toString()
    {
        return "HlaSequenceLoci(allele=" + mAllele + ", sequence=" + sequence() + ')';
    }

    public final boolean consistentWith(final PhasedEvidence evidence)
    {
        int[] nArray = evidence.getAminoAcidIndices();
        return consistentWithAny(evidence.getEvidence().keySet(), Arrays.copyOf(nArray, nArray.length));
    }

    public final boolean consistentWithAny(final Collection<String> targetSequence, final int[] targetIndices)
    {
        // TODO
        return true;
    }

    public final boolean consistentWithAny(final Collection<String> targetSequence, final List<Integer> targetIndices)
    {
        // TODO
        return true;
        
        /*
        boolean bl;
        block3:
        {
            Iterable $receiver$iv = targetSequence;
            if(((Collection) $receiver$iv).isEmpty())
            {
                bl = false;
            }
            else
            {
                for(Object element$iv : $receiver$iv)
                {
                    String it = (String) element$iv;
                    boolean bl2 = false;
                    if(!consistentWith(it, Arrays.copyOf(targetIndices, targetIndices.length)))
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
        
         */
    }

    public final boolean consistentWith(final String targetSequence, final Set<Integer> targetIndices)
    {
        return match(targetSequence, targetIndices) != HlaSequenceMatch.NONE;
    }

    // TODO - require both methods?
    public final boolean consistentWith(final String targetSequence, final int[] targetIndices)
    {
        return false; // match(targetSequence, targetIndices) != HlaSequenceMatch.NONE;
    }

    public final HlaSequenceMatch match(final String targetSequence, final Set<Integer> targetIndices)
    {
        // TODO
        /*
        //int[] nArray = targetIndices;

        if(targetIndices.isEmpty())
            return HlaSequenceMatch.NONE;

        String hlaSequence = sequence(targetIndices);

        if(hlaSequence.length() != targetSequence.length())
            return HlaSequenceMatch.NONE;

        int wildCardCount = 0;
        int n = 0;
        int n2 = ((CharSequence) targetSequence).length();
        while(n < n2)
        {
            void index;
            char target = targetSequence.charAt((int) index);
            if(hlaSequence.charAt((int) index) != '*' && hlaSequence.charAt((int) index) != target)
            {
                return HlaSequenceMatch.NONE;
            }
            if(hlaSequence.charAt((int) index) == '*')
            {
                ++wildCardCount;
            }
            ++index;
        }
        if(wildCardCount > 0)
        {
            if(wildCardCount == targetIndices.length)
            {
                return HlaSequenceMatch.WILD;
            }
            return HlaSequenceMatch.PARTIAL;
        }

         */
        return HlaSequenceMatch.FULL;
    }

    public final List<HlaSequenceLoci> create(final List<HlaSequence> sequences)
    {
        // CHECK
        final String reference = sequences.get(0).getRawSequence();
        return sequences.stream().map(x -> create(mAllele, x.getRawSequence(), reference)).collect(Collectors.toList());
    }

    public static HlaSequenceLoci create(final HlaAllele allele, final String sequence, final String reference)
    {
        final List<String> sequences = Lists.newArrayList();

        int insLength = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            char seqChar = sequence.charAt(i);
            boolean isBaseIgnored = (seqChar == '.' && reference.charAt(i) == '.') || (seqChar == '|');
            boolean isBaseInserted = seqChar != '.' && (i >= reference.length() || reference.charAt(i) == '.');

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

            return new HlaSequenceLoci(allele, sequences);
        }

        // TODO
        return null;
    }
}
