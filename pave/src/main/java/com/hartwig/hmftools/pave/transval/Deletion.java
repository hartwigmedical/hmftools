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

import org.jetbrains.annotations.NotNull;

class Deletion extends ProteinVariant
{
    public Deletion(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @VisibleForTesting
    List<ChangeContext> fittingChanges(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        AminoAcidSequence targetSequence = changeContext.applyDeletion();
        int maxMoves = changeContext.StartPositionInExon;
        if (Transcript.negStrand())
        {
            maxMoves = changeContext.ContainingExon.inExonLength() - changeContext.FinishPositionInExon - 1;
        }
        maxMoves = Math.min(maxMoves, 32);
        List<ChangeContext> result = new ArrayList<>();
        for (int i=0; i<=maxMoves; i++)
        {
            int excisionStart = Transcript.posStrand() ? changeContext.StartPositionInExon - i : changeContext.StartPositionInExon + i;
            int excisionEnd = excisionStart + RefLength * 3 - 1;
            ChangeContext change = new ChangeContext(changeContext.ContainingExon, excisionStart, excisionEnd, Transcript.posStrand(), 0);
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
//        ChangeContext leftMostChange = Transcript.posStrand() ? changes.get(changes.size() - 1) : changes.get(0); CATGTATGAAGTACAGTGGAA
        ChangeContext leftMostChange = changes.get(changes.size() - 1);

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
