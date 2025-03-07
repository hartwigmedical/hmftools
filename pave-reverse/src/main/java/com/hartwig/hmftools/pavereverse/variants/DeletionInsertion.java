package com.hartwig.hmftools.pavereverse.variants;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.SplitCodonSequence;

import org.jetbrains.annotations.NotNull;

public class DeletionInsertion extends ProteinVariant
{
    @NotNull
    public final AminoAcidSequence mAlt;

    public DeletionInsertion(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange,
            @NotNull final AminoAcidSequence alt)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        this.mAlt = alt;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDeletedSection = positionOfFirstAlteredCodon();
        int endOfDeletedSection = startOfDeletedSection + this.mRefLength - 1;
        return completeReferenceAminoAcidSequence().replaceRange(startOfDeletedSection, endOfDeletedSection, mAlt);
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        SplitCodonSequence seq = context.basesForProteinChange(positionOfFirstAlteredCodon(), mRefLength);
        Set<String> destinationBasesForwardStrand = new NucleotidesCalculator(mAlt,seq.retainedPrefix(), seq.retainedSuffix()).allPossibleBaseSequences();
        if(mGene.reverseStrand())
        {
            destinationBasesForwardStrand = destinationBasesForwardStrand.stream().map(Nucleotides::reverseComplementBases).collect(Collectors.toSet());
        }
        int changeStart = seq.locationOfDeletedBases();
        String referenceBases = mGene.forwardStrand()
                ? seq.segmentThatIsModified()
                : Nucleotides.reverseComplementBases(seq.segmentThatIsModified());
        Set<ChangeResult> result = new HashSet<>();
        destinationBasesForwardStrand.forEach(bases -> {
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

    @Override
    boolean seekResultsInCompanionContext(boolean resultsFoundAlready)
    {
        return !resultsFoundAlready;
    }
}
