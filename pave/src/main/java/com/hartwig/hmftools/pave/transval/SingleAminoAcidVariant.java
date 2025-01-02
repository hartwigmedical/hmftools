package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;

public class SingleAminoAcidVariant
{
    @NotNull
    final GeneData gene;
    @NotNull
    final TranscriptData transcript;
    @NotNull
    final TranscriptAminoAcids aminoAcidSequence;
    int position;
    @NotNull final String variant;

    private static boolean isValidAminoAcidName(String s)
    {
        // TODO do we need to handle Z and B???
        return AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s);
    }

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            @NotNull final String variant)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < transcript.length());
        Preconditions.checkArgument(isValidAminoAcidName(variant));
        this.gene = gene;
        this.transcript = transcript;
        this.aminoAcidSequence = aminoAcidSequence;
        this.position = position;
        this.variant = variant;
    }

    public String referenceAminoAcid()
    {
        return aminoAcidSequence.AminoAcids.substring(position - 1, position);
    }

    public String variantAminoAcid()
    {
        return variant;
    }

    public Set<String> possibleVariantCodons()
    {
        return new HashSet<>(AminoAcids.AMINO_ACID_TO_CODON_MAP.get(variant));
    }

    public String referenceCodon(RefGenomeInterface refGenomeSource)
    {
        int codonPosition = 3 * (position - 1);
        CodonExons exonData = exonsForCodonPosition(codonPosition);
        String nucleotides = exonData.retrieveCodon(gene.Chromosome, refGenomeSource);
        if (transcript.negStrand())
        {
            return Nucleotides.reverseComplementBases(nucleotides);
        }
        return nucleotides;
    }

    List<Integer> exonCodingLengths()
    {
        List<ExonData> exons = transcript.exons();
        ArrayList<Integer> result = new ArrayList<>(exons.size());
        for (int i = 0; i < exons.size(); ++i) {
            ExonData exon = exons.get(i);
            int start = codingStartForExon(i);
            int end = i == exons.size() - 1 ? transcript.CodingEnd : exon.End;
            int codingLength = end - start + 1;
            result.add(codingLength);
        }
        if (transcript.negStrand())
        {
            // todo fix this
            List<Integer> reversed = new ArrayList<>(result);
            Collections.reverse(reversed);
            return reversed;
        }
        return result;
    }

    private Integer codingStartForExon(int exonIndex)
    {
        return exonIndex == 0 ? transcript.CodingStart : transcript.exons().get(exonIndex).Start;
    }

    private CodonExons exonsForCodonPosition(int codonPosition)
    {
        List<ExonData> exonsInOrder = new ArrayList<>(transcript.exons());
        if (transcript.negStrand())
        {
            Collections.reverse(exonsInOrder);
        }
//        printExons(exonsInOrder);
        List<Integer> exonLengths = exonCodingLengths();
        int lengthIncludingCurrent = 0;
        for (int i=0; i< exonLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            ExonData exon = exonsInOrder.get(i);
            lengthIncludingCurrent += exonLengths.get(i);
            if (lengthIncludingCurrent > codonPosition)
            {
                ExonData nextExon = (i < exonsInOrder.size() - 1) ? exonsInOrder.get(i + 1) : null;
                int positionOfCodonInCurrentExon = codingStartForExon(i) + (codonPosition - lengthUpToCurrent);
                return new CodonExons(positionOfCodonInCurrentExon, exon, nextExon);
            }
        }
        throw new IllegalArgumentException("No exon found for codon " + codonPosition);
    }

    private void printExons(List<ExonData> exons)
    {
        exons.forEach(this::printExon);
    }

    private void printExon(ExonData exon)
    {
        System.out.println(exon.Rank + " : " + exon.baseLength() + " : "+ exon.Start + "-" + exon.End);
    }
}

