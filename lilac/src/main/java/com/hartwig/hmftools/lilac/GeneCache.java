package com.hartwig.hmftools.lilac;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene_;

import org.apache.commons.lang3.NotImplementedException;

public class GeneCache
{
    public final Map<HlaGene_, TranscriptData> GeneTranscriptMap_;
    public final List<TranscriptData> Transcripts_;

    // TODO: Load the appropriate genes
    public final List<HlaGene_> GeneNames; // long names matching Ensembl

    public final int ExpectAlleleCount;

    public final Map<HlaGene_, List<Integer>> AminoAcidExonBoundaries;
    public final int MaxCommonAminoAcidExonBoundary;

    public final Map<HlaGene_, List<Integer>> NucleotideExonBoundaries;
    public final Map<HlaGene_, Integer> NucleotideLengths;

    public GeneCache(final Map<HlaGene_, TranscriptData> hlaTranscriptMap_)
    {
        GeneTranscriptMap_ = hlaTranscriptMap_;

        // establish other properties and commonly used constants
        GeneNames = GeneTranscriptMap_.keySet().stream().sorted().toList(); // long names matching Ensembl

        Transcripts_ = Lists.newArrayListWithExpectedSize(hlaTranscriptMap_.size());
        GeneNames.forEach(x -> Transcripts_.add(GeneTranscriptMap_.get(x)));

        ExpectAlleleCount = GeneNames.size() * 2;

        AminoAcidExonBoundaries = Maps.newHashMap();
        NucleotideExonBoundaries = Maps.newHashMap();
        NucleotideLengths = Maps.newHashMap();

        for(HlaGene_ geneName : GeneNames)
            setExonBoundaryValues(geneName, GeneTranscriptMap_.get(geneName));

        MaxCommonAminoAcidExonBoundary = findMaxCommonAminoAcidBoundary();
    }

    private void setExonBoundaryValues(final HlaGene_ geneName, final TranscriptData transcriptData_)
    {
        boolean forwardStrand = transcriptData_.posStrand();

        int exonCount = transcriptData_.exons().size();
        int exonIndex = forwardStrand ? 0 : exonCount - 1;

        int codingStart = transcriptData_.CodingStart;
        int codingEnd = transcriptData_.CodingEnd;

        int codingBaseIndex = -1;
        int totalCodingBases = 0;

        List<Integer> aminoAcidExonBoundaries = Lists.newArrayListWithExpectedSize(10);
        AminoAcidExonBoundaries.put(geneName, aminoAcidExonBoundaries);

        List<Integer> nucleotideExonBoundaries = Lists.newArrayListWithExpectedSize(10);
        NucleotideExonBoundaries.put(geneName, nucleotideExonBoundaries);

        while(exonIndex >= 0 && exonIndex < exonCount)
        {
            ExonData exon_ = transcriptData_.exons().get(exonIndex);

            boolean withinCoding = false;

            if(forwardStrand)
            {
                if(exon_.Start > codingEnd)
                    break;

                if(exon_.End >= codingStart)
                    withinCoding = true;
            }
            else
            {
                if(exon_.End < codingStart)
                    break;

                if(exon_.Start <= codingEnd)
                    withinCoding = true;
            }

            // within the coding region
            if(withinCoding)
            {
                if(codingBaseIndex < 0)
                {
                    codingBaseIndex = 0;
                }

                int codingBasesInExon = min(exon_.End, codingEnd) - max(exon_.Start, codingStart) + 1;

                totalCodingBases += codingBasesInExon;

                boolean hasExonBoundary = forwardStrand ? exon_.End <= codingEnd : exon_.Start >= codingStart;

                if(hasExonBoundary)
                {
                    int aminoAcidIndex = totalCodingBases / 3;
                    aminoAcidExonBoundaries.add(aminoAcidIndex);

                    nucleotideExonBoundaries.add(aminoAcidIndex * 3);
                    nucleotideExonBoundaries.add(aminoAcidIndex * 3 + 1);
                    nucleotideExonBoundaries.add(aminoAcidIndex * 3 + 2);
                }
            }

            exonIndex += forwardStrand ? 1 : -1;
        }

        NucleotideLengths.put(geneName, totalCodingBases);
    }

    private int findMaxCommonAminoAcidBoundary()
    {
        int maxCommonAminoAcidBoundary = 0;

        List<List<Integer>> aminoAcidBoundaries = AminoAcidExonBoundaries.values().stream().toList();

        List<Integer> firstSet = aminoAcidBoundaries.get(0);

        for(Integer aaExonBoundary : firstSet)
        {
            boolean allPresent = true;

            for(int i = 1; i < aminoAcidBoundaries.size(); ++i)
            {
                List<Integer> nextSet = aminoAcidBoundaries.get(i);

                if(!nextSet.contains(aaExonBoundary))
                {
                    allPresent = false;
                    break;
                }
            }

            if(allPresent)
            {
                maxCommonAminoAcidBoundary = aaExonBoundary;
            }
        }

        return maxCommonAminoAcidBoundary;
    }
}
