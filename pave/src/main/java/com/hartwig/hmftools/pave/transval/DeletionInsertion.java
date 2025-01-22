package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class DeletionInsertion extends ProteinVariant
{
    final private int DeletionLength;
    @NotNull
    final private String Alt;

    public DeletionInsertion(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            final int deletionLength,
            @NotNull final String variant)
    {
        super(gene, transcript, aminoAcidSequence, position);
        Preconditions.checkArgument(!variant.isBlank());
        for(char ch : variant.toCharArray())
        {
            Preconditions.checkArgument(isValidAminoAcidName("" + ch));
        }
        this.Alt = variant;
        this.DeletionLength = deletionLength;
    }

    @VisibleForTesting
    public SplitSequence referenceBases(RefGenomeInterface refGenome)
    {
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            lengthIncludingCurrent += regionLengths.get(i);
            ChrBaseRegion exon = CodingRegions.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < CodingRegions.size() - 1) ? CodingRegions.get(i + 1) : null;
                if(Transcript.negStrand())
                {
                    //                        int positionOfCodonInCurrentExon = exon.end() - (codonPosition - lengthUpToCurrent);
                    //                        return new CodonRegions(positionOfCodonInCurrentExon, exon, nextExon, false);
                    return null;
                }
                else
                {
                    int positionOfStartInCurrent = exon.start() + (codonPosition - lengthUpToCurrent);
                    int positionOfEndRelativeToStart = positionOfStartInCurrent + DeletionLength * 3 - 1;
                    if(positionOfEndRelativeToStart <= exon.end())
                    {
                        String partInCurrent =
                                refGenome.getBaseString(Gene.Chromosome, positionOfStartInCurrent, positionOfEndRelativeToStart);
                        return new SplitSequence(partInCurrent, null);
                    }
                    else
                    {
                        if(nextExon == null)
                        {
                            throw new IllegalStateException("Variant does not fit given exons.");
                        }
                        String partInCurrentExon = refGenome.getBaseString(Gene.Chromosome, positionOfStartInCurrent, exon.end());
                        int lengthInNext = positionOfEndRelativeToStart - exon.end();
                        String partInNext = refGenome.getBaseString(Gene.Chromosome, nextExon.start(), nextExon.start() + lengthInNext - 1);
                        return new SplitSequence(partInCurrentExon, partInNext);
                    }
                }
            }
        }
        //            throw new IllegalArgumentException("No exon found for codon " + codonPosition);

        return null;
    }

    @Override
    public TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
//        CodonRegions codonRegions = exonsForCodonPosition(positionOfFirstAlteredCodon());
        //        codonRegions.

        /*
        We only consider delinses within a single exon.
        There could conceivably be up to 3 exons involved:
        ...---|---|-] [--|---|---|-] [--|---|---....
        In this example we could replace the entire middle exon resulting in a change
        to the AA that begins in the first exon and finishes in the middle exon
        and a change to the AA that begins in the middle exon and ends in the last one.
         */
        return new TransvalComplexInsertionDeletion(
                Transcript.TransName,
                Gene.Chromosome,
                0,
                false,
                ""
        );
    }

    public String altAminoAcidSequence()
    {
        return Alt;
    }

    @Override
    int changedReferenceSequenceLength()
    {
        return DeletionLength;
    }
}
