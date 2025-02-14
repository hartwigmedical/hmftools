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
    final GeneData mGene;
    @NotNull
    final TranscriptData mTranscript;
    @NotNull
    final TranscriptAminoAcids mAminoAcidSequence;
    private final int mPositionOfFirstAlteredCodon;
    @NotNull
    final List<ChrBaseRegion> mCodingRegions;
    final int mRefLength;

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
    }

    public String referenceAminoAcids()
    {
        return mAminoAcidSequence.AminoAcids.substring(
                mPositionOfFirstAlteredCodon - 1, mPositionOfFirstAlteredCodon + mRefLength - 1);
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
        ChangeContextBuilder exonBuilder = window.seekExonLocation(regionLengths);
        return exonBuilder.build(mGene.Chromosome, genome, mCodingRegions, mTranscript.posStrand());
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

    @VisibleForTesting
    Collection<ChangeResult> findLeftmostApplicableChanges(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        ChangeResult firstExampleResult = applyChange(changeContext).iterator().next();
        int maxMoves = changeContext.StartPositionInExon;
        if(mTranscript.negStrand())
        {
            maxMoves = changeContext.mExon.inExonLength() - changeContext.FinishPositionInExon - 1;
        }
        maxMoves = Math.min(maxMoves, 32);
        Map<String, ChangeResult> results = new HashMap<>();
        int numberOfMovesWithoutAFittingVariant = 0;
        for(int i = 0; i <= maxMoves; i++)
        {
            int start = mTranscript.posStrand() ? changeContext.StartPositionInExon - i : changeContext.StartPositionInExon + i;
            int end = start + mRefLength * 3 - 1;
            ChangeContext change = new ChangeContext(changeContext.mExon, start, end, mTranscript.posStrand(), 0);
            Set<ChangeResult> changesAtThisPosition = applyChange(change);
            boolean fittingVariantFound = false;
            for(ChangeResult changeAtThisPosition : changesAtThisPosition)
            {
                if(firstExampleResult.mAminoAcids.equals(changeAtThisPosition.mAminoAcids))
                {
                    fittingVariantFound = true;
                    results.put(changeAtThisPosition.mBases, changeAtThisPosition);
                    if (singleValueOnlyRequiredForEachStep())
                    {
                        break;
                    }
                }
            }
            if(fittingVariantFound)
            {
                numberOfMovesWithoutAFittingVariant = 0;
            }
            else
            {
//                numberOfMovesWithoutAFittingVariant += 1;
                break;
            }
            if (numberOfMovesWithoutAFittingVariant > 5)
            {
                break;
            }
        }
        return results.values();
    }

    @VisibleForTesting
    boolean singleValueOnlyRequiredForEachStep()
    {
        return false;
    }

    TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        Collection<ChangeResult> changes = findLeftmostApplicableChanges(refGenome);
//        AminoAcidSequence requiredSequence = variantSequence();
//        if(requiredSequence != null)
//        {
//            System.out.println(requiredSequence.sequence());
//            for(ChangeContext change : changes)
//            {
//                ChangeResult cr = change.applyDuplication();
//                System.out.println(cr.mAminoAcids);
//                Preconditions.checkArgument(requiredSequence.sequence().contains(cr.mAminoAcids.sequence()));
//                System.out.println("------- ok: " + cr.mAminoAcids);
//            }
//        }
        if(changes.isEmpty())
        {
            return null;
        }
        Set<TransvalHotspot> hotspots = new HashSet<>();
        changes.forEach(change -> hotspots.add(convertToHotspot(change)));
        return new TransvalVariant(
                mTranscript,
                mGene.Chromosome,
                false,
                hotspots);
    }

    TransvalHotspot convertToHotspot(ChangeResult changeResult)
    {
        return new TransvalHotspot(changeResult.mRefBases, changeResult.altBases, mGene.Chromosome, changeResult.mLocation);
    }

    CodonRegions exonsForCodonPosition(int codonPosition)
    {
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            ChrBaseRegion exon = mCodingRegions.get(i);
            lengthIncludingCurrent += regionLengths.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < mCodingRegions.size() - 1) ? mCodingRegions.get(i + 1) : null;
                if(mTranscript.negStrand())
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
