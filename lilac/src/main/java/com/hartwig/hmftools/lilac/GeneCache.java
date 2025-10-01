package com.hartwig.hmftools.lilac;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public class GeneCache
{
    public final Map<HlaGene, TranscriptData> GeneTranscriptMap;
    public final List<TranscriptData> Transcripts;

    public final List<HlaGene> GeneNames; // long names matching Ensembl

    public final int ExpectAlleleCount;

    public final Map<HlaGene, List<Integer>> AminoAcidExonBoundaries;
    public final int MaxCommonAminoAcidExonBoundary;

    public final Map<HlaGene, List<Integer>> NucleotideExonBoundaries;
    public final Map<HlaGene, Integer> NucleotideLengths;

    public GeneCache(final Map<HlaGene, TranscriptData> hlaTranscriptMap)
    {
        GeneTranscriptMap = hlaTranscriptMap;

        // establish other properties and commonly used constants
        GeneNames = GeneTranscriptMap.keySet().stream().sorted().toList(); // long names matching Ensembl

        Transcripts = Lists.newArrayListWithExpectedSize(hlaTranscriptMap.size());
        GeneNames.forEach(x -> Transcripts.add(GeneTranscriptMap.get(x)));

        ExpectAlleleCount = GeneNames.size() * 2;

        AminoAcidExonBoundaries = Maps.newHashMap();
        NucleotideExonBoundaries = Maps.newHashMap();
        NucleotideLengths = Maps.newHashMap();

        for(HlaGene geneName : GeneNames)
            setExonBoundaryValues(geneName, GeneTranscriptMap.get(geneName));

        MaxCommonAminoAcidExonBoundary = findMaxCommonAminoAcidBoundary();
    }

    private void setExonBoundaryValues(final HlaGene geneName, final TranscriptData transcriptData)
    {
        boolean forwardStrand = transcriptData.posStrand();

        int exonCount = transcriptData.exons().size();
        int exonIndex = forwardStrand ? 0 : exonCount - 1;

        int codingStart = transcriptData.CodingStart;
        int codingEnd = transcriptData.CodingEnd;

        int codingBaseIndex = -1;
        int totalCodingBases = 0;

        List<Integer> aminoAcidExonBoundaries = Lists.newArrayListWithExpectedSize(10);
        AminoAcidExonBoundaries.put(geneName, aminoAcidExonBoundaries);

        List<Integer> nucleotideExonBoundaries = Lists.newArrayListWithExpectedSize(10);
        NucleotideExonBoundaries.put(geneName, nucleotideExonBoundaries);

        while(exonIndex >= 0 && exonIndex < exonCount)
        {
            ExonData exon = transcriptData.exons().get(exonIndex);

            boolean withinCoding = false;

            if(forwardStrand)
            {
                if(exon.Start > codingEnd)
                    break;

                if(exon.End >= codingStart)
                    withinCoding = true;
            }
            else
            {
                if(exon.End < codingStart)
                    break;

                if(exon.Start <= codingEnd)
                    withinCoding = true;
            }

            // within the coding region
            if(withinCoding)
            {
                if(codingBaseIndex < 0)
                {
                    codingBaseIndex = 0;
                }

                int codingBasesInExon = min(exon.End, codingEnd) - max(exon.Start, codingStart) + 1;

                totalCodingBases += codingBasesInExon;

                boolean hasExonBoundary = forwardStrand ? exon.End <= codingEnd : exon.Start >= codingStart;

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
