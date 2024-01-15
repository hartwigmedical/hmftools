package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacConstants.NUC_LENGTH_A;
import static com.hartwig.hmftools.lilac.LilacConstants.NUC_LENGTH_B;
import static com.hartwig.hmftools.lilac.LilacConstants.NUC_LENGTH_C;
import static com.hartwig.hmftools.lilac.LilacUtils.arrayToList;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.read.ReadRecord;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.LociPosition;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class NucleotideFragmentFactory
{
    private final LinkedHashMap<HlaSequenceLoci,SuffixTree> mInsertSuffixTrees;
    private final LinkedHashMap<HlaSequenceLoci,SuffixTree> mDeleteSuffixTrees;
    private final int mMinBaseQuality;
    private final LociPosition mLociPosition;

    public NucleotideFragmentFactory(
            int minBaseQuality, final List<HlaSequenceLoci> inserts, final List<HlaSequenceLoci> deletes, final LociPosition lociPosition)
    {
        mMinBaseQuality = minBaseQuality;
        mLociPosition = lociPosition;

        mInsertSuffixTrees = Maps.newLinkedHashMap();
        mDeleteSuffixTrees = Maps.newLinkedHashMap();

        inserts.stream().forEach(x -> mInsertSuffixTrees.put(x, new SuffixTree(x.sequence())));
        deletes.stream().forEach(x -> mDeleteSuffixTrees.put(x, new SuffixTree(x.sequence())));
    }

    public final Fragment createFragment(final ReadRecord record, final String geneName, final byte geneStrand)
    {
        if(record.ReadStart < 0 || record.ReadEnd < record.ReadStart)
            return null;

        boolean reverseStrand = geneStrand == NEG_STRAND;

        int samCodingStartLoci = reverseStrand
                ? mLociPosition.calcNucelotideLocus(record.PositionEnd) : mLociPosition.calcNucelotideLocus(record.PositionStart);

        int samCodingEndLoci = reverseStrand
                ? mLociPosition.calcNucelotideLocus(record.PositionStart) : mLociPosition.calcNucelotideLocus(record.PositionEnd);

        final char[] codingRegionRead = record.codingRegionRead(reverseStrand);
        final int[] codingRegionQuality = record.codingRegionQuality(reverseStrand);

        if(record.containsIndel() || record.containsSoftClip())
        {
            List<Integer> aminoAcidIndices = calcAminoAcidIndices(samCodingStartLoci, samCodingEndLoci);
            int firstAAIndex = aminoAcidIndices.get(0);
            int nucleotideStartLoci = firstAAIndex * 3;
            String sequence = String.valueOf(codingRegionRead);
            int startLoci = nucleotideStartLoci - samCodingStartLoci;

            if(startLoci < 0 || startLoci >= sequence.length())
            {
                // likely due to a delete in this region
                LL_LOGGER.trace("invalid startLoci({}) requested: read({}) gene({})", startLoci, record.Id, geneName);
                return null;
            }
            String aminoAcids = Codons.aminoAcidFromBases(sequence.substring(startLoci));

            if(!aminoAcids.isEmpty())
            {
                int matchRangeAllowedStart = firstAAIndex - record.SoftClippedStart / 3 - record.maxIndelSize();
                int matchRangeAllowedEnd = firstAAIndex + record.maxIndelSize() + record.SoftClippedEnd / 3;

                Fragment matchedFragment = checkMatchedInsertDeleteSequence(
                        record, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mInsertSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;

                matchedFragment = checkMatchedInsertDeleteSequence(
                        record, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mDeleteSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;
            }

            if(record.containsIndel())
                return null;
        }

        if(samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        List<Integer> lociRange = formRange(samCodingStartLoci, samCodingEndLoci);
        List<String> nucleotides = arrayToList(codingRegionRead);
        List<Integer> qualities = arrayToList(codingRegionQuality);

        return new Fragment(record, geneName, Sets.newHashSet(geneName), lociRange, qualities, nucleotides);
    }

    private Fragment checkMatchedInsertDeleteSequence(
            final ReadRecord record, final String geneName,
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
            if(!fragment.getNucleotideLoci().isEmpty())
                return fragment;
        }

        return null;
    }

    private Fragment createIndelFragment(
            final ReadRecord record, final String geneName, final int startLoci,
            final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        int endLoci = endLoci(startLoci, bamSequence, hlaSequence);
        List<Integer> aminoAcidLoci = formRange(startLoci, endLoci);
        List<Integer> nucleotideLoci = expandIndices(aminoAcidLoci);

        List<String> nucleotides = Lists.newArrayList();

        aminoAcidLoci.stream()
                .map(x -> hlaSequence.sequence(x))
                .map(x -> createNucleotidesFromAminoAcid(x))
                .forEach(x -> nucleotides.addAll(x));

        List<Integer> qualities = nucleotideLoci.stream().map(x -> mMinBaseQuality).collect(Collectors.toList());

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

    public Fragment createAlignmentFragments(final ReadRecord record, final String geneName, final byte geneStrand)
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

    public int calculatePercentileBaseQuality(final List<Fragment> fragments, double percentile)
    {
        int maxBaseQual = mMinBaseQuality * 2; // for purpose of data capture only
        int[] baseQualFrequeny = new int[maxBaseQual + 1];
        long totalBases = 0;

        for(Fragment fragment : fragments)
        {
            for(Integer baseQual : fragment.getRawNucleotideQuality())
            {
                ++totalBases;
                ++baseQualFrequeny[min(baseQual, maxBaseQual)];
            }
        }

        // calculate median
        long percentileEntry = (long)round(totalBases * percentile);
        long cumulativeTotal = 0;

        for(int i = 0; i < baseQualFrequeny.length; ++i)
        {
            cumulativeTotal += baseQualFrequeny[i];

            if(cumulativeTotal >= percentileEntry)
                return i;
        }

        return mMinBaseQuality;
    }

    public static Map<String,int[]> calculateGeneCoverage(final List<Fragment> fragments)
    {
        final Map<String,int[]> geneBaseDepth = Maps.newHashMap();

        geneBaseDepth.put(HLA_A, new int[NUC_LENGTH_A]);
        geneBaseDepth.put(HLA_B, new int[NUC_LENGTH_B]);
        geneBaseDepth.put(HLA_C, new int[NUC_LENGTH_C]);

        for(Fragment fragment : fragments)
        {
            int[] baseDepth = geneBaseDepth.get(fragment.readGene());

            for(int locus : fragment.getRawNucleotideLoci())
            {
                if(locus < baseDepth.length)
                    ++baseDepth[locus];
            }
        }

        return geneBaseDepth;
    }

}
