package com.hartwig.hmftools.pavereverse.protein;

import static java.util.stream.Collectors.toSet;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.SplitCodonSequence;
import com.hartwig.hmftools.pavereverse.util.PRUtils;

public class DeletionInsertion extends ProteinVariant
{
    public final AminoAcidSequence Alt;

    public DeletionInsertion(
            GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange,
            AminoAcidSequence alt)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        Alt = alt;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDeletedSection = positionOfFirstAlteredCodon();
        int endOfDeletedSection = startOfDeletedSection + this.RefLength - 1;
        return completeReferenceAminoAcidSequence().replaceRange(startOfDeletedSection, endOfDeletedSection, Alt);
    }

    @Override
    int numberOfLeftShiftsToTry(final ChangeContext changeContext)
    {
        return 2;
    }

    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        SplitCodonSequence seq = context.basesForProteinChange(positionOfFirstAlteredCodon(), RefLength);
        Set<String> newBases = possibleInsertedNucleotideSequences(seq.retainedPrefix(), seq.retainedSuffix());
        if(reverseStrand())
        {
            newBases = newBases.stream().map(Nucleotides::reverseComplementBases).collect(toSet());
        }
        int changeStart = seq.locationOfDeletedBases();
        String referenceBases = forwardStrand()
                ? seq.segmentThatIsModified()
                : Nucleotides.reverseComplementBases(seq.segmentThatIsModified());
        Set<ChangeResult> result = new HashSet<>();
        newBases.forEach(bases ->
        {
            DeletionInsertionChange changeCanonicaliser = new DeletionInsertionChange(referenceBases, bases);
            String ref = changeCanonicaliser.deleted();
            String alt = changeCanonicaliser.inserted();
            int positionAdjusted = changeStart + changeCanonicaliser.positionOfDeletion();
            String exonBasesAfterChange = context.exonBasesWithReplacementAppliedAtStrandLocation(positionAdjusted, ref, alt);
            AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(exonBasesAfterChange);
            result.add(new ChangeResult(acids, exonBasesAfterChange, positionAdjusted, ref, alt));
        });

        return result;
    }

    @VisibleForTesting
    Set<String> possibleInsertedNucleotideSequences(String prefix, String suffix)
    {
        final NucleotidesCalculator nucleotidesCalculator = new NucleotidesCalculator(Alt, prefix, suffix);
        if(Alt.length() > 2)
        {
            return Set.of(nucleotidesCalculator.anyBaseSequence());
        }
        return nucleotidesCalculator.allPossibleBaseSequences();
    }

    @Override
    boolean seekResultsInCompanionContext(boolean resultsFoundAlready)
    {
        return !resultsFoundAlready;
    }

    @Override
    Collection<ChangeResult> selectChangesToReport(Collection<ChangeResult> changes)
    {
        if(RefLength == 2 && Alt.length() == 2)
        {
            final int[] minEditDistanceSoFar = { Integer.MAX_VALUE };
            Set<ChangeResult> result = new HashSet<>();
            changes.forEach(change ->
            {
                int editDistance = PRUtils.substitutionDistance(change.RefBases, change.AltBases);
                if(editDistance < minEditDistanceSoFar[0])
                {
                    result.clear();
                    result.add(change);
                    minEditDistanceSoFar[0] = editDistance;
                }
                else if(editDistance == minEditDistanceSoFar[0])
                {
                    result.add(change);
                }
            });
            return result;
        }
        return changes;
    }
}
