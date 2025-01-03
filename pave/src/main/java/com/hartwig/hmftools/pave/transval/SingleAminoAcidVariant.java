package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class SingleAminoAcidVariant
{
    @NotNull
    final GeneData Gene;
    @NotNull
    final TranscriptData Transcript;
    @NotNull
    final TranscriptAminoAcids AminoAcidSequence;
    int Position;
    @NotNull
    final String Variant;
    @NotNull
    final List<ChrBaseRegion> CodingRegions;

    private static boolean isValidAminoAcidName(String s)
    {
        // TODO do we need to handle Z and B???
        return AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s);
    }

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            @NotNull final String variant)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < transcript.length());
        Preconditions.checkArgument(isValidAminoAcidName(variant));
        this.Gene = gene;
        this.Transcript = transcript;
        this.AminoAcidSequence = aminoAcidSequence;
        this.Position = position;
        this.Variant = variant;
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

    public String referenceAminoAcid()
    {
        return AminoAcidSequence.AminoAcids.substring(Position - 1, Position);
    }

    public String variantAminoAcid()
    {
        return Variant;
    }

    public Set<String> possibleVariantCodons()
    {
        return new HashSet<>(AminoAcids.AMINO_ACID_TO_CODON_MAP.get(Variant));
    }

    public String referenceCodon(RefGenomeInterface refGenomeSource)
    {
        int codonPosition = 3 * (Position - 1);
        CodonRegions exonData = exonsForCodonPosition(codonPosition);
        return exonData.retrieveCodon(refGenomeSource);
    }

    List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
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
