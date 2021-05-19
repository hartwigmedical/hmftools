package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacUtils.arrayToList;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragment.expandIndices;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragment.merge;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.read.SAMCodingRecord;
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

    public final NucleotideFragment createFragment(final SAMCodingRecord record, final NamedBed codingRegion)
    {
        boolean reverseStrand = record.ReverseStrand;

        int samCodingStartLoci = reverseStrand
                ? mLociPosition.nucelotideLoci(record.PositionEnd) : mLociPosition.nucelotideLoci(record.PositionStart);

        int samCodingEndLoci = reverseStrand
                ? mLociPosition.nucelotideLoci(record.PositionStart) : mLociPosition.nucelotideLoci(record.PositionEnd);

        final char[] codingRegionRead = record.codingRegionRead(reverseStrand);
        final int[] codingRegionQuality = record.codingRegionQuality(reverseStrand);

        if(record.containsIndel() || record.containsSoftClip())
        {
            List<Integer> aminoAcidIndices = calcAminoAcidIndices(samCodingStartLoci, samCodingEndLoci);
            int firstAAIndex = aminoAcidIndices.get(0);
            int nucleotideStartLoci = firstAAIndex * 3;
            String sequence = String.valueOf(codingRegionRead);
            String aminoAcids = Codons.aminoAcids(sequence.substring(nucleotideStartLoci - samCodingStartLoci));

            if(!aminoAcids.isEmpty())
            {
                int matchRangeAllowedStart = firstAAIndex - record.SoftClippedStart / 3 - record.maxIndelSize();
                int matchRangeAllowedEnd = firstAAIndex + record.maxIndelSize() + record.SoftClippedEnd / 3;

                NucleotideFragment matchedFragment = checkMatchedInsertDeleteSequence(
                        record, codingRegion, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mInsertSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;

                matchedFragment = checkMatchedInsertDeleteSequence(
                        record, codingRegion, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mDeleteSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;
            }

            if(record.containsIndel())
                return null;
        }

        if (samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        List<Integer> lociRange = formRange(samCodingStartLoci, samCodingEndLoci);
        List<String> nucleotides = arrayToList(codingRegionRead);
        List<Integer> qualities = arrayToList(codingRegionQuality);

        return new NucleotideFragment(record.Id, Sets.newHashSet(codingRegion.name()), lociRange, qualities, nucleotides);
    }

    private NucleotideFragment checkMatchedInsertDeleteSequence(
            final SAMCodingRecord record, final NamedBed codingRegion,
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

        if(!matchedIndicesList.isEmpty())
        {
            HlaSequenceLoci seqLoci = matchedSeqLoci.get(0);
            List<Integer> filteredAaIndices = matchedIndicesList.get(0);
            NucleotideFragment fragment = createNucleotideSequence(record.Id, codingRegion, filteredAaIndices.get(0), aminoAcids, seqLoci);
            return !fragment.getNucleotideLoci().isEmpty() ? fragment : null;
        }

        return null;
    }

    private NucleotideFragment createNucleotideSequence(
            final String id, final NamedBed codingRegion, final int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
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

        return new NucleotideFragment(id, Sets.newHashSet(codingRegion.name()), nucleotideLoci, qualities, nucleotides);
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

        String codons = Codons.codons(aminoAcid);
        return Lists.newArrayList(
                String.valueOf(codons.charAt(0)), String.valueOf(codons.charAt(1)), codons.substring(2));
    }

    public NucleotideFragment createAlignmentFragments(final SAMCodingRecord record, final NamedBed codingRegion)
    {
        List<NucleotideFragment> fragments = record.alignmentsOnly().stream()
                .filter(x -> x != null)
                .map(x -> createFragment(record, codingRegion)).collect(Collectors.toList());

        if(fragments.isEmpty())
            return null;

        if(fragments.size() == 1)
            return fragments.get(0);

        return NucleotideFragment.reduceById(fragments).get(0);
    }

}
