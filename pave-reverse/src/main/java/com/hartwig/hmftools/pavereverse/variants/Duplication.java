package com.hartwig.hmftools.pavereverse.variants;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.PaddedExon;

public class Duplication extends ProteinVariant
{
    public Duplication(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @Override
    Set<ChangeResult> applyChange(ChangeContext changeContext)
    {
        int changeStart = changeContext.StartPositionInExon;
        int changeEnd = changeContext.FinishPositionInExon + 1;
        PaddedExon exon = changeContext.Exon;
        String bases = exon.baseSequenceWithDuplicationApplied(changeStart, changeEnd, changeContext.IsPositiveStrand);
        AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(bases);
        String duplicated = changeContext.refBases();
        String refBase = duplicated.substring(0, 1);
        return Set.of(new ChangeResult(acids, bases, changeContext.positionOfChangeStartInStrand(), refBase, duplicated));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDuplicatedSection = positionOfFirstAlteredCodon() - 1;
        int endOfDuplicatedSection = startOfDuplicatedSection + this.RefLength;
        return completeReferenceAminoAcidSequence().duplicateRange(startOfDuplicatedSection, endOfDuplicatedSection);
    }
}
