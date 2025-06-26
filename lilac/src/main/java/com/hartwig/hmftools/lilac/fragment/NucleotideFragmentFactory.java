package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacUtils.calcNucelotideLocus;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.utils.Nucleotide;

public class NucleotideFragmentFactory
{
    private final ReferenceData mReferenceData;
    private final LinkedHashMap<HlaSequenceLoci,SuffixTree> mInsertSuffixTrees;
    private final LinkedHashMap<HlaSequenceLoci,SuffixTree> mDeleteSuffixTrees;
    private final byte mMinBaseQuality;

    public NucleotideFragmentFactory(final ReferenceData referenceData)
    {
        mMinBaseQuality = LOW_BASE_QUAL_THRESHOLD;
        mReferenceData = referenceData;

        mInsertSuffixTrees = Maps.newLinkedHashMap();
        mDeleteSuffixTrees = Maps.newLinkedHashMap();

        mReferenceData.AminoAcidSequencesWithInserts.stream().forEach(x -> mInsertSuffixTrees.put(x, new SuffixTree(x.sequence())));
        mReferenceData.AminoAcidSequencesWithDeletes.stream().forEach(x -> mDeleteSuffixTrees.put(x, new SuffixTree(x.sequence())));
    }

    public final Fragment createFragment(final Read read, final String geneName, final byte geneStrand)
    {
        if(read.ReadIndexStart < 0 || read.ReadIndexEnd < read.ReadIndexStart || read.PositionEnd < read.PositionStart)
            return null;

        boolean reverseStrand = geneStrand == NEG_STRAND;

        int codingPositionStartLoci = calcNucelotideLocus(GENE_CACHE.Transcripts, read.PositionStart);
        int codingPositionEndLoci = calcNucelotideLocus(GENE_CACHE.Transcripts, read.PositionEnd);

        int samCodingStartLoci = !reverseStrand ? codingPositionStartLoci : codingPositionEndLoci;
        int samCodingEndLoci = !reverseStrand ? codingPositionEndLoci : codingPositionStartLoci;

        int readLength = read.ReadIndexEnd - read.ReadIndexStart + 1;
        final char[] codingRegionReadBases = new char[readLength];
        final byte[] codingRegionQualities = new byte[readLength];

        read.populateCodingRegion(codingRegionReadBases, codingRegionQualities, reverseStrand);

        if(read.containsIndel() || read.containsSoftClip())
        {
            List<Integer> aminoAcidIndices = calcAminoAcidIndices(samCodingStartLoci, samCodingEndLoci);
            int firstAAIndex = aminoAcidIndices.get(0);
            int nucleotideStartLoci = firstAAIndex * 3;
            String sequence = String.valueOf(codingRegionReadBases);
            int startLoci = nucleotideStartLoci - samCodingStartLoci;

            if(startLoci < 0 || startLoci >= sequence.length())
            {
                // likely due to a delete in this region
                LL_LOGGER.trace("invalid startLoci({}) requested: read({}) gene({})", startLoci, read.Id, geneName);
                return null;
            }
            String aminoAcids = Codons.aminoAcidFromBases(sequence.substring(startLoci));

            if(!aminoAcids.isEmpty())
            {
                int matchRangeAllowedStart = firstAAIndex - read.SoftClippedStart / 3 - read.maxIndelSize();
                int matchRangeAllowedEnd = firstAAIndex + read.maxIndelSize() + read.SoftClippedEnd / 3;

                Fragment matchedFragment = checkMatchedInsertDeleteSequence(
                        read, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mInsertSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;

                matchedFragment = checkMatchedInsertDeleteSequence(
                        read, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mDeleteSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;
            }

            if(read.containsIndel())
                return null;
        }

        if(samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        int rangeLength = samCodingEndLoci - samCodingStartLoci + 1;
        List<Nucleotide> nucleotides = Lists.newArrayListWithCapacity(rangeLength);

        for(int i = 0; i < rangeLength; ++i)
        {
            nucleotides.add(new Nucleotide(
                    samCodingStartLoci + i, Byte.valueOf(codingRegionQualities[i]), String.valueOf(codingRegionReadBases[i])));
        }

        return new Fragment(read, geneName, Sets.newHashSet(geneName), nucleotides);
    }

    private Fragment checkMatchedInsertDeleteSequence(
            final Read record, final String geneName,
            final String aminoAcids, int matchRangeAllowedStart, int matchRangeAllowedEnd,
            final LinkedHashMap<HlaSequenceLoci,SuffixTree> sequenceMap)
    {
        List<List<Integer>> matchedIndicesList = Lists.newArrayList();
        List<HlaSequenceLoci> matchedSeqLoci = Lists.newArrayList();

        for(Map.Entry<HlaSequenceLoci,SuffixTree> entry : sequenceMap.entrySet())
        {
            HlaSequenceLoci seqLoci = entry.getKey();
            List<Integer> aaIndices = entry.getValue().indices(aminoAcids);

            List<Integer> filteredAaIndices = aaIndices.stream()
                    .filter(x -> x >= matchRangeAllowedStart && x <= matchRangeAllowedEnd)
                    .collect(Collectors.toList());

            if(!filteredAaIndices.isEmpty())
            {
                matchedSeqLoci.add(seqLoci);
                matchedIndicesList.add(filteredAaIndices);
            }
        }

        for(int i = 0; i < matchedSeqLoci.size(); ++i)
        {
            HlaSequenceLoci seqLoci = matchedSeqLoci.get(i);
            List<Integer> filteredAaIndices = matchedIndicesList.get(i);
            Fragment fragment = createIndelFragment(record, geneName, filteredAaIndices.get(0), aminoAcids, seqLoci);
            if(!fragment.nucleotidesByLoci().isEmpty())
                return fragment;
        }

        return null;
    }

    private Fragment createIndelFragment(
            final Read record, final String geneName, final int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        int endLoci = endLoci(startLoci, bamSequence, hlaSequence);
        List<Integer> aminoAcidLoci = formRange(startLoci, endLoci);
        List<Integer> nucleotideLoci = expandIndices(aminoAcidLoci);

        List<String> nucleotides = Lists.newArrayList();

        aminoAcidLoci.stream()
                .map(x -> hlaSequence.sequence(x))
                .map(x -> createNucleotidesFromAminoAcid(x))
                .forEach(x -> nucleotides.addAll(x));

        List<Byte> qualities = nucleotideLoci.stream().map(x -> mMinBaseQuality).collect(Collectors.toList());

        return new Fragment(record, geneName, Sets.newHashSet(geneName), nucleotideLoci, qualities, nucleotides);
    }

    private int endLoci(int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        String tmpSequence = "";

        for(int locus = startLoci; locus < hlaSequence.length(); ++locus)
        {
            tmpSequence += hlaSequence.getSequences().get(locus);

            if(!bamSequence.startsWith(tmpSequence))
                return locus - 1;
        }

        return hlaSequence.length() - 1;
    }

    public static List<String> createNucleotidesFromAminoAcid(final String aminoAcid)
    {
        if(aminoAcid.equals(DEL_STR))
            return Lists.newArrayList(DEL_STR, DEL_STR, DEL_STR);

        String codons = Codons.aminoAcidsToCodons(aminoAcid);
        return Lists.newArrayList(
                String.valueOf(codons.charAt(0)), String.valueOf(codons.charAt(1)), codons.substring(2));
    }

    public Fragment createAlignmentFragments(final Read record, final String geneName, final byte geneStrand)
    {
        List<Fragment> fragments = record.alignmentsOnly().stream()
                .filter(x -> x != null)
                .map(x -> createFragment(record, geneName, geneStrand))
                .filter(x -> x != null)
                .collect(Collectors.toList());

        if(fragments.isEmpty())
            return null;

        if(fragments.size() == 1)
            return fragments.get(0);

        return FragmentUtils.mergeFragmentsById(fragments).get(0);
    }

    public byte calculatePercentileBaseQuality(final List<Fragment> fragments, double percentile)
    {
        // calculates the nth percentile base quality for each fragment's nucleotides
        int maxBaseQual = mMinBaseQuality * 2;
        int[] baseQualFrequeny = new int[maxBaseQual + 1];
        long totalBases = 0;

        for(Fragment fragment : fragments)
        {
            for(byte baseQual : Nucleotide.qualities(fragment.rawNucleotidesByLoci().values()))
            {
                ++totalBases;
                ++baseQualFrequeny[min(baseQual, maxBaseQual)];
            }
        }

        long percentileEntry = round(totalBases * percentile);
        long cumulativeTotal = 0;

        for(int i = 0; i < baseQualFrequeny.length; ++i)
        {
            cumulativeTotal += baseQualFrequeny[i];

            if(cumulativeTotal >= percentileEntry)
                return (byte) i;
        }

        return mMinBaseQuality;
    }

    public static Map<String,int[]> calculateGeneCoverage(final List<Fragment> fragments)
    {
        final Map<String,int[]> geneBaseDepth = Maps.newHashMap();

        for(String geneName : GENE_CACHE.GeneNames)
        {
            int geneNucleotideCount = GENE_CACHE.NucleotideLengths.get(geneName);
            geneBaseDepth.put(geneName, new int[geneNucleotideCount]);
        }

        for(Fragment fragment : fragments)
        {
            int[] baseDepth = geneBaseDepth.get(fragment.readGene());

            for(int locus : fragment.rawNucleotidesByLoci().keySet())
            {
                if(locus < baseDepth.length)
                    ++baseDepth[locus];
            }
        }

        return geneBaseDepth;
    }

}
