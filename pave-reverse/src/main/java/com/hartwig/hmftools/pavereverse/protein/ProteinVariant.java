package com.hartwig.hmftools.pavereverse.protein;

import static java.util.stream.Collectors.toSet;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.ChangeContextBuilder;

public abstract class ProteinVariant extends HgvsVariant
{
    final TranscriptAminoAcids AminoAcidsTranscript;
    private final int mPositionOfFirstAlteredCodon;
    public final int RefLength;
    private final AminoAcidSequence mCompleteReferenceSequence;

    ProteinVariant(
            GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            final int positionOfFirstAlteredCodon,
            final int refSequenceLength)
    {
        super(gene, transcript);
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(positionOfFirstAlteredCodon >= 0);
        Preconditions.checkArgument(positionOfFirstAlteredCodon < transcript.length());
        AminoAcidsTranscript = aminoAcidSequence;
        mPositionOfFirstAlteredCodon = positionOfFirstAlteredCodon;
        RefLength = refSequenceLength;
        mCompleteReferenceSequence = AminoAcidSequence.parse(AminoAcidsTranscript.AminoAcids);
    }

    AminoAcidSequence completeReferenceAminoAcidSequence()
    {
        return mCompleteReferenceSequence;
    }

    public int positionOfFirstAlteredCodon()
    {
        return mPositionOfFirstAlteredCodon;
    }

    ChangeContext getChangeContext(RefGenomeInterface genome)
    {
        CodonWindow window = new CodonWindow(positionOfFirstAlteredCodon(), RefLength);
        List<Integer> regionLengths = codingRegionLengths();
        ChangeContextBuilder changeContextBuilder = window.seekExonLocation(regionLengths);
        return changeContextBuilder.build(chromosome(), genome, codingRegions(), forwardStrand());
    }

    abstract Set<ChangeResult> applyChange(ChangeContext changeContext);

    abstract AminoAcidSequence variantSequence();

    @VisibleForTesting
    public final AminoAcidSequence replaceExonAminoAcids(int exon, AminoAcidSequence replacement)
    {
        if(exon == 0 && (replacement.length() == 0))
        {
            return AminoAcidSequence.empty();
        }
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
        if(reverseStrand())
        {
            result = changeContext.Exon.inExonLength() - changeContext.FinishPositionInExon - 1;
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
        if(changeContext.companionContext() != null && seekResultsInCompanionContext(!result.isEmpty()))
        {
            result.addAll(findLeftmostApplicableChangesForContext(changeContext.companionContext()));
        }
        return result;
    }

    boolean seekResultsInCompanionContext(boolean resultsFoundAlready)
    {
        return true;
    }

    boolean isConsistentWithThisVariant(AminoAcidSequence aminoAcidSequence)
    {
        return aminoAcidSequence.equals(variantSequence());
    }

    private Collection<ChangeResult> findLeftmostApplicableChangesForContext(ChangeContext changeContext)
    {
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
                AminoAcidSequence changedAminoAcidSequence = replaceExonAminoAcids(changeContext.Exon.Index, changeAtThisPosition.Acids);
                if(isConsistentWithThisVariant(changedAminoAcidSequence))
                {
                    results.put(changeAtThisPosition.ExonBases, changeAtThisPosition);
                }
            }
        }
        return new HashSet<>(results.values());
    }

    public BaseSequenceVariants calculateVariant(RefGenomeInterface refGenome)
    {
        Collection<ChangeResult> changes = findLeftmostApplicableChanges(refGenome);
        Collection<ChangeResult> changesToReport = selectChangesToReport(changes);
        Set<BaseSequenceChange> hotspots = changesToReport.stream().map(cr -> cr.asChange(chromosome(), refGenome)).collect(toSet());
        return new BaseSequenceVariants(transcriptData(), chromosome(), hotspots);
    }

    Collection<ChangeResult> selectChangesToReport(Collection<ChangeResult> changes)
    {
        return changes;
    }

    boolean doesNotMatchVariantUpToLastAminoAcid(AminoAcidSequence candidate)
    {
        AminoAcidSequence variant = variantSequence();
        if(candidate.length() < variant.length() + 1)
        {
            return true;
        }
        for(int i = 0; i < variant.length(); i++)
        {
            if(!candidate.get(i).equals(variant.get(i)))
            {
                return true;
            }
        }
        return false;
    }
}
