package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacUtils.arrayToList;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragment.merge;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.read.SAMCodingRecord;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.LociPosition;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class NucleotideFragmentFactory
{
    private final Map<HlaSequenceLoci,SuffixTree> mInsertSuffixTrees;
    private final Map<HlaSequenceLoci,SuffixTree> mDeleteSuffixTrees;
    private final int mMinBaseQuality;
    private final LociPosition mLociPosition;

    public NucleotideFragmentFactory(
            int minBaseQuality, final List<HlaSequenceLoci> inserts, final List<HlaSequenceLoci> deletes, final LociPosition lociPosition)
    {
        mMinBaseQuality = minBaseQuality;
        mLociPosition = lociPosition;

        mInsertSuffixTrees = Maps.newHashMap();
        mDeleteSuffixTrees = Maps.newHashMap();

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

                // CHECK !
                for(Map.Entry<HlaSequenceLoci,SuffixTree> entry : mInsertSuffixTrees.entrySet())
                {
                    List<Integer> aaIndices = entry.getValue().indices(aminoAcids);
                    List<Integer> filteredAaIndices = aaIndices.stream()
                            .filter(x -> x >= matchRangeAllowedStart && x <= matchRangeAllowedEnd)
                            .collect(Collectors.toList());

                    if(!filteredAaIndices.isEmpty())
                    {
                        return createNucleotideSequence(record.Id, codingRegion, filteredAaIndices.get(0), aminoAcids, entry.getKey());
                    }
                }

                // TODO - turn these 2 into function
                for(Map.Entry<HlaSequenceLoci,SuffixTree> entry : mDeleteSuffixTrees.entrySet())
                {
                    List<Integer> aaIndices = entry.getValue().indices(aminoAcids);
                    List<Integer> filteredAaIndices = aaIndices.stream()
                            .filter(x -> x >= matchRangeAllowedStart && x <= matchRangeAllowedEnd)
                            .collect(Collectors.toList());

                    if(!filteredAaIndices.isEmpty())
                    {
                        return createNucleotideSequence(record.Id, codingRegion, filteredAaIndices.get(0), aminoAcids, entry.getKey());
                    }
                }
            }

            if(record.containsIndel())
                return null;
        }

        if (samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        // CHECK
        List<Integer> lociRange = formRange(samCodingStartLoci, samCodingEndLoci);
        List<String> nucleotides = arrayToList(codingRegionRead);
        List<Integer> qualities = arrayToList(codingRegionQuality);

        return new NucleotideFragment(record.Id, Sets.newHashSet(codingRegion.name()), lociRange, qualities, nucleotides);
    }

    private final NucleotideFragment createNucleotideSequence(
            final String id, final NamedBed codingRegion, final int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        int endLoci = endLoci(startLoci, bamSequence, hlaSequence);
        List<Integer> aminoAcidLoci = formRange(startLoci, endLoci);
        List<Integer> nucleotideLoci = Lists.newArrayList();
        aminoAcidLoci.stream().forEach(x -> nucleotideLoci.addAll(Lists.newArrayList(3 * x, 3 * x + 1, 3 * x + 2)));

        List<String> nucleotides = Lists.newArrayList();

        aminoAcidLoci.stream()
                .map(x -> hlaSequence.sequence(x))
                .map(x -> createNucleotidesFromAminoAcid(x))
                .forEach(x -> nucleotides.addAll(x));

        List<Integer> qualities = nucleotideLoci.stream().map(x -> mMinBaseQuality).collect(Collectors.toList());

        return new NucleotideFragment(id, Sets.newHashSet(codingRegion.name()), nucleotideLoci, qualities, nucleotides);
    }

    private final int endLoci(int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        String tmpSequence = "";

        // CHECK end of for-loop
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
        if(aminoAcid.equals("."))
            return Lists.newArrayList(".", ".", ".");

        String codons = Codons.codons(aminoAcid);
        return Lists.newArrayList(
                String.valueOf(codons.charAt(0)), String.valueOf(codons.charAt(0)), codons.substring(2));
    }

    public NucleotideFragment createAlignmentFragments(final SAMCodingRecord record, final NamedBed codingRegion)
    {
        List<NucleotideFragment> fragments = record.alignmentsOnly().stream()
                .filter(x -> x != null)
                .map(x -> createFragment(record, codingRegion)).collect(
                        Collectors.toList());

        // CHECK: usage in SAMReader - implies only 2?!?
        if(fragments.isEmpty() || fragments.size() != 2)
            return null;

        return merge(fragments.get(0), fragments.get(1));
    }

}
