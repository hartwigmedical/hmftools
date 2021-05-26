package com.hartwig.hmftools.lilac.seq;

import static com.hartwig.hmftools.lilac.LilacConstants.getExonGeneBoundries;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILDCARD;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DELETION;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.EXON_BOUNDARY;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.IDENTICAL;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class HlaSequenceLoci
{
    public final HlaAllele Allele;

    private final List<String> mSequences;

    public HlaSequenceLoci(final HlaAllele allele, final List<String> sequences)
    {
        Allele = allele;
        mSequences = Lists.newArrayList();
        mSequences.addAll(sequences);
    }

    public int length() { return mSequences.size(); }

    public List<String> getSequences() { return mSequences; }

    public boolean containsInserts() { return mSequences.stream().anyMatch(x -> x.length() > 1); }

    public boolean containsDeletes() { return mSequences.stream().anyMatch(x -> x.equals(DEL_STR)); }

    public boolean containsIndels() { return containsDeletes() || containsInserts(); }

    public boolean containsWildcards() { return mSequences.stream().anyMatch(x -> x.equals(WILD_STR)); }

    public String sequence(int locus) { return mSequences.get(locus); }

    public String sequence(int startLocus, int endLocus)
    {
        StringJoiner sj = new StringJoiner("");
        for(int i = startLocus; i <= endLocus; ++i)
        {
            sj.add(mSequences.get(i));
        }

        return sj.toString();
    }

    public String sequence(final List<Integer> indices)
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

    public String sequence()
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
        return String.format("HlaSequenceLoci(allele=%s, sequence=%s)", Allele, sequence());
    }

    public boolean consistentWith(final PhasedEvidence evidence)
    {
        return consistentWithAny(
                evidence.getEvidence().keySet().stream().collect(Collectors.toList()), evidence.getAminoAcidIndices());
    }

    public boolean consistentWithAny(final Set<String> targetSequences, final int[] targetIndices)
    {
        return consistentWithAny(targetSequences.stream().collect(Collectors.toList()), targetIndices);
    }

    public boolean consistentWithAny(final List<String> targetSequences, final int[] targetIndices)
    {
        List<Integer> tiList = Lists.newArrayList();
        Arrays.stream(targetIndices).forEach(x -> tiList.add(x));
        return consistentWithAny(targetSequences, tiList);
    }

    public boolean consistentWithAny(final List<String> targetSequences, final List<Integer> targetIndices)
    {
        return targetSequences.stream().anyMatch(x -> determineMatchType(x, targetIndices) != SequenceMatchType.NONE);
    }

    public boolean consistentWith(final String targetSequence, final int[] targetIndices)
    {
        List<Integer> tiList = Lists.newArrayList();
        Arrays.stream(targetIndices).forEach(x -> tiList.add(x));
        return determineMatchType(targetSequence, tiList) != SequenceMatchType.NONE;
    }

    public static List<Integer> filterExonBoundaryWildcards(final HlaSequenceLoci sequence, final List<Integer> loci)
    {
        // exclude any wildcard location matching an exon boundary
        if(!sequence.containsWildcards())
            return loci;

        List<Integer> aminoAcidExonBoundaries = getExonGeneBoundries(sequence.Allele.Gene);

        return loci.stream()
                .filter(x -> !aminoAcidExonBoundaries.contains(x) || !sequence.sequence(x).equals(WILD_STR))
                .collect(Collectors.toList());
    }

    public static List<Integer> filterWildcards(final HlaSequenceLoci sequence, final List<Integer> loci)
    {
        // exclude any wildcard location matching an exon boundary
        if(!sequence.containsWildcards())
            return loci;

        return loci.stream()
                .filter(x -> !sequence.sequence(x).equals(WILD_STR))
                .collect(Collectors.toList());
    }

    public SequenceMatchType determineMatchType(final String targetSequence, final List<Integer> targetIndices)
    {
        if(targetIndices.isEmpty())
            return SequenceMatchType.NONE;

        String hlaSequence = sequence(targetIndices);

        if(hlaSequence.length() != targetSequence.length())
            return SequenceMatchType.NONE;

        int wildCardCount = 0;

        for(int i = 0; i < targetSequence.length(); ++i)
        {
            char targetChr = targetSequence.charAt(i);
            char hlaSeqChr = hlaSequence.charAt(i);

            if(hlaSeqChr != WILDCARD && hlaSeqChr != targetChr)
                return SequenceMatchType.NONE;

            if (hlaSeqChr == WILDCARD)
                wildCardCount++;
        }

        if (wildCardCount > 0)
            return SequenceMatchType.WILD;

        return SequenceMatchType.FULL;
    }

    public static HlaSequenceLoci create(final HlaAllele allele, final String sequence, final String reference)
    {
        final List<String> sequences = Lists.newArrayList();

        int insLength = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            char seqChar = sequence.charAt(i);
            char refChar = i < reference.length() ? reference.charAt(i) : IDENTICAL;
            boolean isBaseInserted = seqChar != DELETION && (i >= reference.length() || refChar == DELETION);
            boolean isBaseIgnored = (seqChar == '.' && i < reference.length() && refChar == DELETION) || (seqChar == EXON_BOUNDARY);

            if(insLength > 0 && !isBaseInserted)
            {
                String insertStr = sequence.substring(i - insLength, i);
                sequences.set(sequences.size() - 1, sequences.get(sequences.size() - 1) + insertStr);
                insLength = 0;
            }

            if(isBaseInserted)
            {
                insLength++;
            }
            else if(!isBaseIgnored)
            {
                // String locusSequence = seqChar; // sequence[i].toString()

                if(seqChar == IDENTICAL)
                    sequences.add(String.valueOf(refChar));
                else
                    sequences.add(String.valueOf(seqChar));
            }
        }

        return new HlaSequenceLoci(allele, sequences);
    }
}
