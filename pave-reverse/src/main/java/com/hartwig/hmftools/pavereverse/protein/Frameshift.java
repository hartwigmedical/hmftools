package com.hartwig.hmftools.pavereverse.protein;

import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcid;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;

import org.apache.commons.lang3.tuple.Pair;

public class Frameshift extends ProteinVariant
{
    public final AminoAcid FirstChangedAminoAcid;

    public Frameshift(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        Preconditions.checkArgument(refRange.length() == 1);
        FirstChangedAminoAcid = refRange.aminoAcidAtStart();
    }

    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        Pair<String, String> baseToLeftAndBaseDeleted = context.forwardStrandBaseAndLeftNeighbour();
        String newBases = context.Exon.baseSequenceWithSingleBaseRemoved(context.StartPositionInExon, context.IsPositiveStrand);
        AminoAcidSequence newAminoAcids = AminoAcidSequence.fromNucleotides(newBases);
        String alt = baseToLeftAndBaseDeleted.getLeft();
        String ref = alt + baseToLeftAndBaseDeleted.getRight();
        int changePosition = context.insertionPoint();
        if(!context.IsPositiveStrand)
        {
            changePosition -= 1;
        }
        return Set.of(new ChangeResult(newAminoAcids, newBases, changePosition, ref, alt));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().subsequenceUpToInclusive(positionOfFirstAlteredCodon() - 1);
    }

    @Override
    boolean isConsistentWithThisVariant(AminoAcidSequence candidate)
    {
        // The sequences need to match up to the change point and then differ
        // at least at that position.
        if(doesNotMatchVariantUpToLastAminoAcid(candidate))
        {
            return false;
        }
        return !candidate.get(variantSequence().length()).equals(FirstChangedAminoAcid);
    }
}
