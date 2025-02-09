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
    final int RefLength;

    public ProteinVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int positionOfFirstAlteredCodon,
            final int refSequenceLength)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(positionOfFirstAlteredCodon >= 0);
        Preconditions.checkArgument(positionOfFirstAlteredCodon < transcript.length());
        this.Gene = gene;
        this.Transcript = transcript;
        this.AminoAcidSequence = aminoAcidSequence;
        this.PositionOfFirstAlteredCodon = positionOfFirstAlteredCodon;
        this.RefLength = refSequenceLength;
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
                PositionOfFirstAlteredCodon - 1, PositionOfFirstAlteredCodon +  RefLength - 1);
    }

    public int positionOfFirstAlteredCodon()
    {
        return PositionOfFirstAlteredCodon;
    }

    public int positionOfLastAlteredCodon()
    {
        return PositionOfFirstAlteredCodon + RefLength - 1;
    }

    @VisibleForTesting
    List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }

    @VisibleForTesting
    ChangeContext getChangeContext2(RefGenomeInterface genome)
    {
        CodonWindow window = new CodonWindow(positionOfFirstAlteredCodon(), RefLength);
        List<Integer> regionLengths = codingRegionLengths();
        ChangeContextBuilder exonBuilder = window.seekExonLocation(regionLengths);
        return exonBuilder.build(Gene.Chromosome, genome, CodingRegions);
    }

    @VisibleForTesting
    ChangeContext getChangeContext(RefGenomeInterface genome)
    {
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        ChrBaseRegion previousExon;
        ChrBaseRegion exon = null;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            int aminoAcidsStartingInPreviousExons = lengthUpToCurrent / 3;
            lengthIncludingCurrent += regionLengths.get(i);
            previousExon = exon;
            int lengthOfFirstCodonOfCurrentExonInPreviousExon = lengthUpToCurrent % 3;
            String basesOfFirstCodonOfCurrentExonInPreviousExon = "";
            exon = CodingRegions.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < CodingRegions.size() - 1) ? CodingRegions.get(i + 1) : null;
                int lengthOfLastCodonOfCurrentExonInCurrentExon = lengthIncludingCurrent % 3;
                int lengthOfLastCodonOfCurrentExonInNextExon = (3 - lengthOfLastCodonOfCurrentExonInCurrentExon) % 3;
                String basesOfLastCodonOfCurrentExonInNextExon = "";
                if(Transcript.negStrand())
                {
                    if (lengthOfFirstCodonOfCurrentExonInPreviousExon > 0) {
                        int start = previousExon.start();
                        int stop = previousExon.start() + lengthOfFirstCodonOfCurrentExonInPreviousExon - 1;
                        basesOfFirstCodonOfCurrentExonInPreviousExon = genome.getBaseString(Gene.Chromosome, start, stop);
                    }
                    if (lengthOfLastCodonOfCurrentExonInNextExon > 0) {
                        int stop = nextExon.end() ;
                        int start = stop - lengthOfLastCodonOfCurrentExonInNextExon + 1;
                        basesOfLastCodonOfCurrentExonInNextExon = genome.getBaseString(Gene.Chromosome, start, stop );
                    }
                    final int relativePositionOfStart = codonPosition - lengthUpToCurrent;
                    int relativePositionOfEnd = relativePositionOfStart + RefLength * 3 - 1;
                    int absolutePositionOfEnd = exon.end() - relativePositionOfEnd;
                    if(absolutePositionOfEnd >= exon.start())
                    {
                        String exonBases = genome.getBaseString(Gene.Chromosome, exon.start(), exon.end());
                        String prefix = genome.getBaseString(Gene.Chromosome, exon.start() - 5, exon.start() - 1);
                        PaddedExon containingExon = new PaddedExon(basesOfFirstCodonOfCurrentExonInPreviousExon, basesOfLastCodonOfCurrentExonInNextExon, exonBases, exon.start(), prefix);
                        return new ChangeContext(containingExon, relativePositionOfStart, relativePositionOfEnd, false, aminoAcidsStartingInPreviousExons + 1);
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    final int relativePositionOfStart = codonPosition - lengthUpToCurrent;
                    final int relativePositionOfEnd = relativePositionOfStart + RefLength * 3 - 1;
                    int absolutePositionOfEnd = exon.start() + relativePositionOfEnd;
                    if(absolutePositionOfEnd <= exon.end())
                    {
                        if (lengthOfFirstCodonOfCurrentExonInPreviousExon > 0) {
                            int start = previousExon.end() - lengthOfFirstCodonOfCurrentExonInPreviousExon + 1;
                            int stop = previousExon.end();
                            aminoAcidsStartingInPreviousExons += 1;
                            basesOfFirstCodonOfCurrentExonInPreviousExon = genome.getBaseString(Gene.Chromosome, start, stop);
                        }
                        if (lengthOfLastCodonOfCurrentExonInNextExon > 0) {
                            int start = nextExon.start();
                            int stop = start + lengthOfLastCodonOfCurrentExonInNextExon - 1 ;
                            basesOfLastCodonOfCurrentExonInNextExon = genome.getBaseString(Gene.Chromosome, start, stop );
                        }
                        String exonBases = genome.getBaseString(Gene.Chromosome, exon.start(), exon.end());
                        String prefix = genome.getBaseString(Gene.Chromosome, exon.start() - 5, exon.start() - 1);
                        PaddedExon containingExon = new PaddedExon(basesOfFirstCodonOfCurrentExonInPreviousExon, basesOfLastCodonOfCurrentExonInNextExon, exonBases, exon.start(), prefix);
                        return new ChangeContext(containingExon, relativePositionOfStart, relativePositionOfEnd, true, aminoAcidsStartingInPreviousExons + 1);
                    }
                    else
                    {
                        return null;
                    }
                }
            }
        }
        return null;
    }

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
