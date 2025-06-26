package com.hartwig.hmftools.lilac.seq;

import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DELETION;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.EXON_BOUNDARY;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.IDENTICAL;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILDCARD;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.FULL;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.WILD;

import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaSequenceLoci
{
    public final HlaAllele Allele;

    private final List<String> mSequences;
    private final boolean mHasDeletes;
    private final boolean mHasInserts;
    private final boolean mHasWildcards;
    private boolean mHasExonBoundaryWildcards;
    private final int mWildcardCount;

    public HlaSequenceLoci(final HlaAllele allele, final List<String> sequences)
    {
        Allele = allele;
        mSequences = Lists.newArrayList();
        mSequences.addAll(sequences);

        mHasDeletes = mSequences.stream().anyMatch(x -> x.equals(DEL_STR));
        mHasInserts = mSequences.stream().anyMatch(x -> x.length() > 1);

        mWildcardCount = (int)mSequences.stream().filter(x -> x.equals(WILD_STR)).count();
        mHasWildcards = mWildcardCount > 0;

        if(mHasWildcards)
            Allele.setHasWildcard(true);

        mHasExonBoundaryWildcards = false;
    }

    public int length() { return mSequences.size(); }

    public List<String> getSequences() { return mSequences; }

    public boolean hasInserts() { return mHasInserts; }
    public boolean hasDeletes() { return mHasDeletes; }
    public boolean hasIndels() { return hasDeletes() || hasInserts(); }
    public boolean hasWildcards() { return mHasWildcards; }
    public int wildcardCount() { return mWildcardCount; }

    public boolean hasExonBoundaryWildcards() { return mHasExonBoundaryWildcards; }

    public void setExonBoundaryWildcards(final List<Integer> exonBoundaries)
    {
        for(Integer locus : exonBoundaries)
        {
            if(locus >= mSequences.size())
                return;

            if(mSequences.get(locus).equals(WILD_STR))
            {
                mHasExonBoundaryWildcards = true;
                return;
            }
        }
    }

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

    public String sequence(final Collection<Integer> indices)
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
                Lists.newArrayList(evidence.getEvidence().keySet()), evidence.getAminoAcidLoci());
    }

    public boolean consistentWithAny(final List<String> targetSequences, final List<Integer> targetIndices)
    {
        return targetSequences.stream().anyMatch(x -> determineMatchType(x, targetIndices) != SequenceMatchType.MISMATCH);
    }

    public boolean consistentWith(final String targetSequence, final List<Integer> targetIndices)
    {
        return determineMatchType(targetSequence, targetIndices) != SequenceMatchType.MISMATCH;
    }

    public static List<Integer> filterExonBoundaryWildcards(final HlaSequenceLoci sequence, final List<Integer> loci)
    {
        // exclude any wildcard location matching an exon boundary
        if(!sequence.hasWildcards())
            return loci;

        List<Integer> aminoAcidExonBoundaries = getAminoAcidExonBoundaries(sequence.Allele.Gene);

        return loci.stream()
                .filter(x -> !aminoAcidExonBoundaries.contains(x)
                        || x < sequence.length() && !sequence.sequence(x).equals(WILD_STR))
                .collect(Collectors.toList());
    }

    public static List<Integer> filterWildcards(final HlaSequenceLoci sequence, final List<Integer> loci)
    {
        // exclude any wildcard location matching an exon boundary
        if(!sequence.hasWildcards())
            return loci;

        return loci.stream()
                .filter(x -> !sequence.sequence(x).equals(WILD_STR))
                .collect(Collectors.toList());
    }

    public SequenceMatchType determineMatchType(final List<String> targetSequences, final Collection<Integer> targetLoci)
    {
        if(targetLoci.isEmpty())
            return SequenceMatchType.MISMATCH;

        boolean hasWildcardMatch = false;

        int i = 0;
        for(Integer locus : targetLoci)
        {
            if(locus >= mSequences.size())
                return SequenceMatchType.MISMATCH;

            String sequenceStr = mSequences.get(locus);

            if(sequenceStr.equals(WILD_STR))
                hasWildcardMatch = true;
            else if(!sequenceStr.equals(targetSequences.get(i)))
                return SequenceMatchType.MISMATCH;

            i++;
        }

        return hasWildcardMatch ? WILD : FULL;
    }

    public SequenceMatchType determineMatchType(final String targetSequence, final Collection<Integer> targetIndices)
    {
        if(targetIndices.isEmpty())
            return SequenceMatchType.MISMATCH;

        String hlaSequence = sequence(targetIndices);

        if(hlaSequence.equals(targetSequence))
            return FULL;

        if(hlaSequence.length() != targetSequence.length())
            return SequenceMatchType.MISMATCH;

        int wildCardCount = 0;

        for(int i = 0; i < targetSequence.length(); ++i)
        {
            char targetChr = targetSequence.charAt(i);
            char hlaSeqChr = hlaSequence.charAt(i);

            if(hlaSeqChr != WILDCARD && hlaSeqChr != targetChr)
                return SequenceMatchType.MISMATCH;

            if (hlaSeqChr == WILDCARD)
                wildCardCount++;
        }

        if (wildCardCount > 0)
            return WILD;

        return FULL;
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

    public static HlaSequenceLoci buildAminoAcidSequenceFromNucleotides(
            final HlaSequenceLoci nucSequence, final HlaSequenceLoci sequenceTemplate)
    {
        StringBuilder sequence = new StringBuilder();
        for(int i = 0; i < nucSequence.length() - 2; i = i + 3)
        {
            String first = nucSequence.getSequences().get(i);
            String second = nucSequence.getSequences().get(i + 1);
            String third = nucSequence.getSequences().get(i + 2);

            if(first.equals(DEL_STR) || second.equals(DEL_STR) || third.equals(DEL_STR))
                sequence.append(DEL_STR);
            else
                sequence.append(Codons.aminoAcidFromBases(first + second + third));
        }

        return create(nucSequence.Allele.asFourDigit(), sequence.toString(), sequenceTemplate.sequence());
    }
}
