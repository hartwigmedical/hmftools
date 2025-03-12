package com.hartwig.hmftools.lilac;

import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;

public class GeneCache
{
    public final Map<String,TranscriptData> GeneTranscriptMap;
    public final List<String> GeneIds; // strips off 'HLA-' prefix, used for logging and look-ups
    public final List<String> GeneNames; // long names matching Ensembl

    public final int ExpectAlleleCount;

    public final List<String> ExcludedAlleles;
    public final List<String> StopLossAlleles;

    public final Map<String,List<Integer>> AminoAcidExonBoundaries;
    public final int MaxCommonAminoAcidExonBoundary;

    public final Map<String,List<Integer>> NucleotideExonBoundaries;
    public final Map<String,Integer> NucleotideLengths;

    /*
    public static final List<String> GENE_IDS = Lists.newArrayList(GENE_A, GENE_B, GENE_C);
    public static final List<String> HLA_GENES = Lists.newArrayList(HLA_A, HLA_B, HLA_C);

    public static final int EXPECTED_ALLELE_COUNT = 6;

    public static final List<String> EXCLUDED_ALLELES = Lists.newArrayList("A*31:135", "A*33:191", "A*02:783", "B*07:282");

    // common INDEL associated with allele C*04:09N
    public static final String STOP_LOSS_ON_C_ALLELE = "C*04:09N";

    public static final List<Integer> A_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 364);
    public static final List<Integer> B_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348);
    public static final List<Integer> C_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 338, 349, 365);

    public static final int NUC_LENGTH_A = 1098;
    public static final int NUC_LENGTH_B = 1089;
    public static final int NUC_LENGTH_C = 1101;

    public static final int MAX_AMINO_ACID_BOUNDARY = 298;

    public static final Map<String,List<Integer>> NUCLEOTIDE_EXON_BOUNDARIES = Maps.newHashMap();

    // common routines using constants
    public static List<Integer> getAminoAcidExonBoundaries(final String gene)
    {
        return gene.equals(GENE_A) ? A_EXON_BOUNDARIES : (gene.equals(GENE_B) ? B_EXON_BOUNDARIES : C_EXON_BOUNDARIES);
    }

    public static List<Integer> getNucleotideExonBoundaries(final String gene)
    {
        return NUCLEOTIDE_EXON_BOUNDARIES.get(gene);
    }

    */

    public GeneCache(
            final Map<String,TranscriptData> hlaTranscriptMap, final List<String> excludedAlleles, final List<String> stopLossAlleles)
    {
        GeneTranscriptMap = hlaTranscriptMap;

        // establish other properties and commonly used constants
        GeneNames = GeneTranscriptMap.keySet().stream().toList(); // long names matching Ensembl
        GeneIds = GeneNames.stream().map(x -> shortGeneName(x)).collect(Collectors.toList());

        ExpectAlleleCount = GeneIds.size() * 2;

        ExcludedAlleles = excludedAlleles;
        StopLossAlleles = stopLossAlleles;

        AminoAcidExonBoundaries = Maps.newHashMap();
        NucleotideExonBoundaries = Maps.newHashMap();
        NucleotideLengths = Maps.newHashMap();

        for(String geneName : GeneNames)
        {
            setExonBoundaryValues(geneName, GeneTranscriptMap.get(geneName));
        }

        MaxCommonAminoAcidExonBoundary = findMaxCommonAminoAcidBoundary();
    }

    private void setExonBoundaryValues(final String geneName, final TranscriptData transcriptData)
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
                    withinCoding = true;;
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
                    nucleotideExonBoundaries.add(totalCodingBases);

                    int aminoAcidIndex = totalCodingBases / 3;
                    aminoAcidExonBoundaries.add(aminoAcidIndex);
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

    public static String shortGeneName(final String gene)
    {
        return gene.startsWith(HLA_PREFIX) ? gene.substring(gene.length() - 1) : gene;
    }

    public static String longGeneName(final String gene)
    {
        return gene.length() == 1 ? HLA_PREFIX + gene : gene;
    }

}
