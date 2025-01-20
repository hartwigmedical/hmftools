package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;


public abstract class ProteinVariant
{
    @NotNull
    final GeneData Gene;
    @NotNull
    final TranscriptData Transcript;
    @NotNull
    final TranscriptAminoAcids AminoAcidSequence;
    int Position;
    @NotNull
    final List<ChrBaseRegion> CodingRegions;

    static boolean isValidAminoAcidName(String s)
    {
        return AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s);
    }


    /*
    SingleAminoAcidVariant
    - position
    - length, always 1
    - alt aa, always length 1

    DeletionInsertion
    - position
    - deletion length
    - inserted values

    Deletion
    - position
    - deletion length

    Insertion
    - position
    - (deletion length = 0)
    - inserted values

    Duplication
    - position
    - duplication length

    Frameshift
    - position
     */

    public ProteinVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < transcript.length());
        this.Gene = gene;
        this.Transcript = transcript;
        this.AminoAcidSequence = aminoAcidSequence;
        this.Position = position;
        List<ChrBaseRegion> codingRegions = transcript.exons().stream()
                .filter(exonData -> exonData.End >= transcript.CodingStart && exonData.Start <= transcript.CodingEnd)
                .map(exonData -> new ChrBaseRegion(gene.Chromosome, Math.max(exonData.Start, transcript.CodingStart), Math.min(exonData.End, transcript.CodingEnd)))
                .collect(Collectors.toList());
        if(Transcript.negStrand())
        {
            List<ChrBaseRegion> reversed = new ArrayList<>(codingRegions);
            Collections.reverse(reversed);
            CodingRegions = Collections.unmodifiableList(reversed);
        }
        else
        {
            CodingRegions = codingRegions;
        }
    }

    public String referenceAminoAcids()
    {
        return AminoAcidSequence.AminoAcids.substring(Position - 1, Position +  changedReferenceSequenceLength() - 1);
    }

    abstract int changedReferenceSequenceLength();
}
