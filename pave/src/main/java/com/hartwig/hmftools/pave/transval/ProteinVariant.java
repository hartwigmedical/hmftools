package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
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

    ChangeContext getChangeContext(RefGenomeInterface genome)
    {
        CodonWindow window = new CodonWindow(positionOfFirstAlteredCodon(), RefLength);
        List<Integer> regionLengths = codingRegionLengths();
        ChangeContextBuilder exonBuilder = window.seekExonLocation(regionLengths);
        return exonBuilder.build(Gene.Chromosome, genome, CodingRegions, Transcript.posStrand());
    }

    ChangeResult applyChange(ChangeContext changeContext)
    {
        return null;
    }

    @VisibleForTesting
    Collection<ChangeContext> findLeftmostApplicableChanges(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        ChangeResult firstExampleResult = applyChange(changeContext);
        int maxMoves = changeContext.StartPositionInExon;
        if (Transcript.negStrand())
        {
            maxMoves = changeContext.ContainingExon.inExonLength() - changeContext.FinishPositionInExon - 1;
        }
        maxMoves = Math.min(maxMoves, 32);
        Map<String,ChangeContext> results = new HashMap<>();
        for (int i=0; i<=maxMoves; i++)
        {
            int excisionStart = Transcript.posStrand() ? changeContext.StartPositionInExon - i : changeContext.StartPositionInExon + i;
            int excisionEnd = excisionStart + RefLength * 3 - 1;
            ChangeContext change = new ChangeContext(changeContext.ContainingExon, excisionStart, excisionEnd, Transcript.posStrand(), 0);
            ChangeResult changeAtThisPosition = applyChange(change);
            if(firstExampleResult.mAminoAcids.equals(changeAtThisPosition.mAminoAcids))
            {
                results.put(changeAtThisPosition.mBases, change);
            }
        }
        return results.values();
    }

    TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        Collection<ChangeContext> changes = findLeftmostApplicableChanges(refGenome);
        if(changes.isEmpty())
        {
            return null;
        }
        Set<TransvalHotspot> hotspots = new HashSet<>();
        changes.forEach(change -> hotspots.add(convertToHotspot(change)));
        return new TransvalVariant(
                Transcript,
                Gene.Chromosome,
                false,
                hotspots);

    }

    abstract TransvalHotspot convertToHotspot(ChangeContext changeContext);

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
