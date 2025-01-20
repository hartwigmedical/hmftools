package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class SingleAminoAcidVariant extends ProteinVariant
{
    @NotNull
    final String Alt;
    @NotNull
    private final CodonRegions RegionsDefiningCodon;

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            @NotNull final String variant)
    {
        super(gene, transcript, aminoAcidSequence, position);
        Preconditions.checkArgument(isValidAminoAcidName(variant));
        this.Position = position;
        this.Alt = variant;
        int codonPosition = 3 * (Position - 1);
        RegionsDefiningCodon = exonsForCodonPosition(codonPosition);
    }

    public String altValue()
    {
        return Alt;
    }

    public String referenceCodon(RefGenomeInterface refGenomeSource)
    {
        return RegionsDefiningCodon.retrieveCodon(refGenomeSource);
    }

    public boolean codonIsInSingleExon()
    {
        return RegionsDefiningCodon.codonIsInSingleExon();
    }

    public CodonRegions regionsDefiningCodon()
    {
        return RegionsDefiningCodon;
    }

    @VisibleForTesting
    List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }

    @Override
    int changedReferenceSequenceLength()
    {
        return 1;
    }

    private CodonRegions exonsForCodonPosition(int codonPosition)
    {
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            ChrBaseRegion exon = CodingRegions.get(i);
            lengthIncludingCurrent += regionLengths.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < CodingRegions.size() - 1) ? CodingRegions.get(i + 1) : null;
                if(Transcript.negStrand())
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
