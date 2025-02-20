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

public class ProteinVariant
{
    @NotNull
    final GeneData mGene;
    @NotNull
    final TranscriptData mTranscript;
    @NotNull
    final TranscriptAminoAcids mAminoAcidSequence;
    private final int mPositionOfFirstAlteredCodon;
    @NotNull
    final List<ChrBaseRegion> mCodingRegions;
    final int mRefLength;
    @NotNull
    private final AminoAcidSequence mCompleteReferenceSequence;

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
        this.mGene = gene;
        this.mTranscript = transcript;
        this.mAminoAcidSequence = aminoAcidSequence;
        this.mPositionOfFirstAlteredCodon = positionOfFirstAlteredCodon;
        this.mRefLength = refSequenceLength;
        List<ChrBaseRegion> codingRegions = transcript.exons().stream()
                .filter(exonData -> exonData.End >= transcript.CodingStart && exonData.Start <= transcript.CodingEnd)
                .map(exonData -> new ChrBaseRegion(gene.Chromosome, Math.max(exonData.Start, transcript.CodingStart), Math.min(exonData.End, transcript.CodingEnd)))
                .collect(Collectors.toList());
        if(mTranscript.negStrand())
        {
            List<ChrBaseRegion> reversed = new ArrayList<>(codingRegions);
            Collections.reverse(reversed);
            mCodingRegions = Collections.unmodifiableList(reversed);
        }
        else
        {
            mCodingRegions = codingRegions;
        }
        mCompleteReferenceSequence = AminoAcidSequence.parse(mAminoAcidSequence.AminoAcids);
    }

    public String referenceAminoAcids()
    {
        return mAminoAcidSequence.AminoAcids.substring(
                mPositionOfFirstAlteredCodon - 1, mPositionOfFirstAlteredCodon + mRefLength - 1);
    }

    AminoAcidSequence referenceAminoAcidSequence()
    {
        return mCompleteReferenceSequence;
    }

    public int positionOfFirstAlteredCodon()
    {
        return mPositionOfFirstAlteredCodon;
    }

    public int positionOfLastAlteredCodon()
    {
        return mPositionOfFirstAlteredCodon + mRefLength - 1;
    }

    @VisibleForTesting
    List<Integer> codingRegionLengths()
    {
        return mCodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }

    ChangeContext getChangeContext(RefGenomeInterface genome)
    {
        CodonWindow window = new CodonWindow(positionOfFirstAlteredCodon(), mRefLength);
        List<Integer> regionLengths = codingRegionLengths();
        ChangeContextBuilder changeContextBuilder = window.seekExonLocation(regionLengths);
        return changeContextBuilder.build(mGene.Chromosome, genome, mCodingRegions, mTranscript.posStrand());
    }

    @NotNull
    Set<ChangeResult> applyChange(ChangeContext changeContext)
    {
        return new HashSet<>();
    }

    AminoAcidSequence variantSequence()
    {
        return null;
    }

    @NotNull
    AminoAcidSequence replaceExonAminoAcids(int exon, @NotNull AminoAcidSequence replacement)
    {
        int lengthUpToExon = 0;
        List<Integer> regionLengths = codingRegionLengths();
        for(int i = 0; i < exon; i++)
        {
            lengthUpToExon += regionLengths.get(i);
        }
        int numberOfAAsBeforeExon = lengthUpToExon / 3;
        int lengthIncludingExon = lengthUpToExon + regionLengths.get(exon);
        int numberOfAAsIncludingExon = lengthIncludingExon / 3;
        if(lengthIncludingExon % 3 != 0)
        {
            numberOfAAsIncludingExon++;
        }
        AminoAcidSequence left = mCompleteReferenceSequence.subsequenceUpToInclusive(numberOfAAsBeforeExon);
        AminoAcidSequence right = mCompleteReferenceSequence.subsequenceAfterExclusive(numberOfAAsIncludingExon - 1);
        return left.append(replacement).append(right);
    }

    int numberOfLeftShiftsToTry(ChangeContext changeContext)
    {
        int result = changeContext.StartPositionInExon;
        if(mTranscript.negStrand())
        {
            result = changeContext.mExon.inExonLength() - changeContext.FinishPositionInExon - 1;
        }
        return Math.min(result, 32);
    }

    @VisibleForTesting
    Collection<ChangeResult> findLeftmostApplicableChanges(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        Set<ChangeResult> result = new HashSet<>(findLeftmostApplicableChangesForContext(changeContext));
        // If the change context has a companion (so it's straddling a splice junction)
        // we need to consider any changes from that context.
        if(changeContext.companionContext() != null)
        {
            result.addAll(findLeftmostApplicableChangesForContext(changeContext.companionContext()));
        }
        return result;
    }

    private Collection<ChangeResult> findLeftmostApplicableChangesForContext(ChangeContext changeContext)
    {
        AminoAcidSequence targetSequence = variantSequence();
        final Set<ChangeResult> changeResults = applyChange(changeContext);
        if(changeResults.isEmpty())
        {
            return changeResults;
        }
        int maxMoves = numberOfLeftShiftsToTry(changeContext);
        Map<String, ChangeResult> results = new HashMap<>();
        for(int i = 0; i <= maxMoves; i++)
        {
            ChangeContext change = changeContext.shiftLeft(i);
            Set<ChangeResult> changesAtThisPosition = applyChange(change);
            for(ChangeResult changeAtThisPosition : changesAtThisPosition)
            {
                AminoAcidSequence effectOfChange = replaceExonAminoAcids(changeContext.mExon.mIndex, changeAtThisPosition.mAminoAcids);
                if(effectOfChange.equals(targetSequence))
                {
                    results.put(changeAtThisPosition.mBases, changeAtThisPosition);
                }
            }
        }
        return new HashSet<>(results.values());
    }

    TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        Collection<ChangeResult> changes = findLeftmostApplicableChanges(refGenome);
        Set<TransvalHotspot> hotspots = new HashSet<>();
        changes.forEach(change -> hotspots.add(change.toHotspot(mGene.Chromosome)));
        return new TransvalVariant(mTranscript, mGene.Chromosome, false, hotspots);
    }
}
