package com.hartwig.hmftools.pavereverse.variants;

import static com.hartwig.hmftools.pavereverse.aa.AminoAcid.START;

import java.util.HashSet;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

public class StartLost extends SingleCodonVariant
{
    public StartLost(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition());
        Preconditions.checkArgument(refRange.length() == 1);
        Preconditions.checkArgument(refRange.startPosition() == 1);
        Preconditions.checkArgument(refRange.aminoAcidAtStart().equals(START));
    }

    @Override
    public AminoAcidSequence variantSequence()
    {
        return AminoAcidSequence.empty();
    }

    @Override
    Set<CodonChange> possibleVariants(CodonWithinExons codon)
    {
        Set<CodonChange> changes = new HashSet<>();
        changes.add(new CodonChange("ATG", "ATA"));
        changes.add(new CodonChange("ATG", "ATC"));
        changes.add(new CodonChange("ATG", "ATT"));
        changes.add(new CodonChange("ATG", "AAG"));
        changes.add(new CodonChange("ATG", "ACG"));
        changes.add(new CodonChange("ATG", "AGG"));
        changes.add(new CodonChange("ATG", "CTG"));
        changes.add(new CodonChange("ATG", "GTG"));
        changes.add(new CodonChange("ATG", "TTG"));
        return changes;
    }

    @Override
    boolean isConsistentWithThisVariant(AminoAcidSequence candidate)
    {
        return candidate.length() == 0;
    }

    @Override
    AminoAcidSequence convertBaseSequence(String bases)
    {
        Preconditions.checkArgument(bases.length() > 3);
        Preconditions.checkArgument(!bases.startsWith("ATG"));
        return AminoAcidSequence.empty();
    }
}
