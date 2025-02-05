package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

class Deletion extends ProteinVariant
{
    final private int DeletionLength;
    @NotNull
    final AminoAcidRange RefRange;

    public Deletion(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition());
        this.RefRange = refRange;
        DeletionLength = refRange.length();
    }

    @VisibleForTesting
    ChangeContext getChangeContext(RefGenomeInterface genome)
    {
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        ChrBaseRegion previousExon = null;
        ChrBaseRegion exon = null;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
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
                        int stop = nextExon.end() - 1;
                        int start = stop - lengthOfLastCodonOfCurrentExonInNextExon;
                        basesOfFirstCodonOfCurrentExonInPreviousExon = genome.getBaseString(Gene.Chromosome, start, stop );
                    }
                    final int relativePositionOfStart = codonPosition - lengthUpToCurrent;
                    int positionOfStartInCurrent = exon.end() - relativePositionOfStart;
                    int relativePositionOfEnd = relativePositionOfStart + DeletionLength * 3 - 1;
                    int absolutePositionOfEnd = exon.end() - relativePositionOfEnd;
                    if(absolutePositionOfEnd >= exon.start())
                    {
                        String exonBases = genome.getBaseString(Gene.Chromosome, exon.start(), exon.end());
//                        exonBases = Nucleotides.reverseComplementBases(exonBases);
                        ExtendedExon containingExon = new ExtendedExon(basesOfFirstCodonOfCurrentExonInPreviousExon, basesOfLastCodonOfCurrentExonInNextExon, exonBases, exon.start());
                        return new ChangeContext(containingExon, relativePositionOfStart, relativePositionOfEnd, false);
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    final int relativePositionOfStart = codonPosition - lengthUpToCurrent;
                    int absolutePositionOfStart = exon.start() + relativePositionOfStart;
                    final int relativePositionOfEnd = relativePositionOfStart + DeletionLength * 3 - 1;
                    int absolutePositionOfEnd = exon.start() + relativePositionOfEnd;
                    if(absolutePositionOfEnd <= exon.end())
                    {
                        String exonBases = genome.getBaseString(Gene.Chromosome, exon.start(), exon.end());
                        ExtendedExon containingExon = new ExtendedExon("", "", exonBases, exon.start());
                        return new ChangeContext(containingExon, relativePositionOfStart, relativePositionOfEnd, true);
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

    @Override
    int changedReferenceSequenceLength()
    {
        return DeletionLength;
    }

    @VisibleForTesting
    List<ChangeContext> fittingChanges(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        AminoAcidSequence targetSequence = changeContext.applyDeletion();
        int deletionStart = changeContext.startPositionInExon;
        List<ChangeContext> result = new ArrayList<>();
        for (int i=0; i<=deletionStart; i++)
        {
            int excisionStart = changeContext.startPositionInExon - i;
            int excisionEnd = excisionStart + DeletionLength * 3 - 1;
            ChangeContext change = new ChangeContext(changeContext.containingExon, excisionStart, excisionEnd, Transcript.posStrand());
            AminoAcidSequence effect = change.applyDeletion();
            if(targetSequence.equals(effect))
            {
                result.add(change);
            }
        }
        return result;
    }

    @Override
    TransvalVariant calculateVariant(final RefGenomeInterface refGenome)
    {
        List<ChangeContext> changes = fittingChanges(refGenome);
        if(changes.isEmpty())
        {
            return null;
        }
        ChangeContext leftMostChange = Transcript.posStrand() ? changes.get(changes.size() - 1) : changes.get(0);
//        ChangeContext leftMostChange = changes.get(changes.size() - 1);

        TransvalHotspot changeHotspot = leftMostChange.hotspot(Gene.Chromosome);
        Set<TransvalHotspot> hotspots = new HashSet<>();
        hotspots.add(leftMostChange.hotspot(Gene.Chromosome));
        return new TransvalVariant(
                Transcript,
                Gene.Chromosome,
                changeHotspot.mPosition,
                false,
                changeHotspot.Ref,
                hotspots);
    }
}
