package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;


public abstract class ProteinVariant
{
    @NotNull
    final GeneData Gene;
    @NotNull
    final TranscriptData Transcript;
    @NotNull
    final TranscriptAminoAcids AminoAcidSequence;
    private final int PositionOfFirstAlteredCodon;
    @NotNull
    final List<ChrBaseRegion> CodingRegions;

    public ProteinVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int positionOfFirstAlteredCodon)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(positionOfFirstAlteredCodon >= 0);
        Preconditions.checkArgument(positionOfFirstAlteredCodon < transcript.length());
        this.Gene = gene;
        this.Transcript = transcript;
        this.AminoAcidSequence = aminoAcidSequence;
        this.PositionOfFirstAlteredCodon = positionOfFirstAlteredCodon;
        List<ChrBaseRegion> codingRegions = transcript.exons().stream()
                .filter(exonData -> exonData.End >= transcript.CodingStart && exonData.Start <= transcript.CodingEnd)
                .map(exonData -> new ChrBaseRegion(gene.Chromosome, Math.max(exonData.Start, transcript.CodingStart), Math.min(exonData.End, transcript.CodingEnd)))
                .collect(Collectors.toList());
        if(Transcript.negStrand())
        {
            List<ChrBaseRegion> reversed = new ArrayList<>(codingRegions);
            Collections.reverse(reversed);
            CodingRegions = Collections.unmodifiableList(reversed);
        }
        else
        {
            CodingRegions = codingRegions;
        }
    }

    public String referenceAminoAcids()
    {
        return AminoAcidSequence.AminoAcids.substring(
                PositionOfFirstAlteredCodon - 1, PositionOfFirstAlteredCodon +  changedReferenceSequenceLength() - 1);
    }

    public int positionOfFirstAlteredCodon()
    {
        return PositionOfFirstAlteredCodon;
    }

    public int positionOfLastAlteredCodon()
    {
        return PositionOfFirstAlteredCodon + changedReferenceSequenceLength() - 1;
    }

    @VisibleForTesting
    List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }

    abstract int changedReferenceSequenceLength();

    abstract TransvalVariant calculateVariant(RefGenomeInterface refGenome);

    CodonRegions exonsForCodonPosition(int codonPosition)
    {
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            ChrBaseRegion exon = CodingRegions.get(i);
            lengthIncludingCurrent += regionLengths.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < CodingRegions.size() - 1) ? CodingRegions.get(i + 1) : null;
                if(Transcript.negStrand())
                {
                    int positionOfCodonInCurrentExon = exon.end() - (codonPosition - lengthUpToCurrent);
                    return new CodonRegions(positionOfCodonInCurrentExon, exon, nextExon, false);
                }
                else
                {
                    int positionOfCodonInCurrentExon = exon.start() + (codonPosition - lengthUpToCurrent);
                    return new CodonRegions(positionOfCodonInCurrentExon, exon, nextExon);
                }
            }
        }
        throw new IllegalArgumentException("No exon found for codon " + codonPosition);
    }
}
