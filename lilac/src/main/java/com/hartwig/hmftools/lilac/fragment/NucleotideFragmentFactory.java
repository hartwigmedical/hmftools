package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.calcNucelotideLocus;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.evidence.Nucleotide.MISSING_BASE_QUAL;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class NucleotideFragmentFactory
{
    private final ReferenceData mReferenceData;
    private final LinkedHashMap<HlaSequenceLoci, SuffixTree> mInsertSuffixTrees;
    private final LinkedHashMap<HlaSequenceLoci, SuffixTree> mDeleteSuffixTrees;

    public NucleotideFragmentFactory(final ReferenceData referenceData)
    {
        mReferenceData = referenceData;

        mInsertSuffixTrees = Maps.newLinkedHashMap();
        mDeleteSuffixTrees = Maps.newLinkedHashMap();

        mReferenceData.AminoAcidSequencesWithInserts.forEach(x -> mInsertSuffixTrees.put(x, new SuffixTree(x.sequence())));
        mReferenceData.AminoAcidSequencesWithDeletes.forEach(x -> mDeleteSuffixTrees.put(x, new SuffixTree(x.sequence())));
    }

    public Fragment createFragment(final Read read, final HlaGene geneName, final byte geneStrand)
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
        final boolean[] codingRegionIsLowQuals = new boolean[readLength];
        final byte[] codingRegionQualities = new byte[readLength];

        read.populateCodingRegion(codingRegionReadBases, codingRegionIsLowQuals, codingRegionQualities, reverseStrand);

        if(read.containsValidIndel() || read.containsSoftClip())
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
                int matchRangeAllowedStart = firstAAIndex - read.SoftClippedStart / 3 - read.maxValidIndelSize();
                int matchRangeAllowedEnd = firstAAIndex + read.maxValidIndelSize() + read.SoftClippedEnd / 3;

                Fragment matchedFragment = checkMatchedInsertDeleteSequence(
                        read, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mInsertSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;

                matchedFragment = checkMatchedInsertDeleteSequence(
                        read, geneName, aminoAcids, matchRangeAllowedStart, matchRangeAllowedEnd, mDeleteSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;
            }

            if(read.containsValidIndel())
                return null;
        }

        if(samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        int rangeLength = samCodingEndLoci - samCodingStartLoci + 1;
        List<Nucleotide> nucleotides = Lists.newArrayListWithCapacity(rangeLength);

        Queue<Indel> indelQueue = new PriorityQueue<>(Comparator.comparingInt((Indel x) -> x.ReadIndex));
        indelQueue.addAll(read.getIgnoredIndels());

        final int inc = reverseStrand ? -1 : 1;
        final int iStart = reverseStrand ? readLength - 1 : 0;
        int samCodingLoci = reverseStrand ? samCodingEndLoci : samCodingStartLoci;
        int readIndex = read.ReadIndexStart;
        for(int i = iStart; i >= 0 && i < readLength; )
        {
            if(!indelQueue.isEmpty() && indelQueue.peek().ReadIndex < readIndex)
                indelQueue.poll();

            Indel currentIndel = null;
            if(!indelQueue.isEmpty() && indelQueue.peek().ReadIndex == readIndex)
                currentIndel = indelQueue.poll();

            if(currentIndel == null)
            {
                nucleotides.add(new Nucleotide(
                        samCodingLoci, codingRegionIsLowQuals[i], codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));

                samCodingLoci += inc;
                readIndex++;
                i += inc;
                continue;
            }

            if(currentIndel.IsInsert)
            {
                nucleotides.add(new Nucleotide(
                        samCodingLoci, codingRegionIsLowQuals[i], codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));

                samCodingLoci += inc;
                readIndex += 1 + currentIndel.Length;
                i += inc * (1 + currentIndel.Length);
                continue;
            }

            // currentIndel is del
            nucleotides.add(new Nucleotide(
                    samCodingLoci, codingRegionIsLowQuals[i], codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));
            samCodingLoci += inc;
            readIndex++;
            i += inc;

            for(int j = 1; j < currentIndel.Ref.length(); j++)
            {
                if(samCodingLoci < samCodingStartLoci || samCodingLoci > samCodingEndLoci)
                    break;

                char base = currentIndel.Ref.charAt(j);
                if(reverseStrand)
                    base = Nucleotides.swapDnaBase(base);

                nucleotides.add(Nucleotide.createHighQual(samCodingLoci, String.valueOf(base)));

                samCodingLoci += inc;
            }
        }

        if(reverseStrand)
            Collections.reverse(nucleotides);

        return new Fragment(read, geneName, Sets.newHashSet(geneName), nucleotides);
    }

    private Fragment checkMatchedInsertDeleteSequence(
            final Read record, final HlaGene geneName,
            final String aminoAcids, int matchRangeAllowedStart, int matchRangeAllowedEnd,
            final LinkedHashMap<HlaSequenceLoci, SuffixTree> sequenceMap)
    {
        List<List<Integer>> matchedIndicesList = Lists.newArrayList();
        List<HlaSequenceLoci> matchedSeqLoci = Lists.newArrayList();

        for(Map.Entry<HlaSequenceLoci, SuffixTree> entry : sequenceMap.entrySet())
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

    private static Fragment createIndelFragment(
            final Read record, final HlaGene geneName, final int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        int endLoci = endLoci(startLoci, bamSequence, hlaSequence);
        List<Integer> aminoAcidLoci = formRange(startLoci, endLoci);
        List<Integer> nucleotideLoci = expandIndices(aminoAcidLoci);

        List<String> nucleotides = Lists.newArrayList();

        aminoAcidLoci.stream()
                .map(hlaSequence::sequence)
                .map(NucleotideFragmentFactory::createNucleotidesFromAminoAcid)
                .forEach(nucleotides::addAll);

        List<Boolean> isLowQuals = nucleotideLoci.stream().map(x -> false).toList();
        return Fragment.createFromIsLowQuals(record, geneName, Sets.newHashSet(geneName), nucleotideLoci, isLowQuals, nucleotides);
    }

    private static int endLoci(int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
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

    public Fragment createAlignmentFragments(final Read record, final HlaGene geneName, final byte geneStrand)
    {
        List<Fragment> fragments = record.alignmentsOnly().stream()
                .filter(Objects::nonNull)
                .map(x -> createFragment(record, geneName, geneStrand))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        if(fragments.isEmpty())
            return null;

        if(fragments.size() == 1)
            return fragments.get(0);

        return FragmentUtils.mergeFragmentsById(fragments).get(0);
    }

    public static byte calculatePercentileBaseQuality(final Iterable<Fragment> fragments, double percentile)
    {
        // calculates the nth percentile base quality for each fragment's nucleotides
        int maxBaseQual = LOW_BASE_QUAL_THRESHOLD * 2;
        int[] baseQualFrequency = new int[maxBaseQual + 1];
        long totalBases = 0;

        for(Fragment fragment : fragments)
        {
            for(byte baseQual : Nucleotide.qualities(fragment.rawNucleotidesByLoci().values()))
            {
                if(baseQual == MISSING_BASE_QUAL)
                    continue;

                ++totalBases;
                ++baseQualFrequency[min(baseQual, maxBaseQual)];
            }
        }

        long percentileEntry = round(totalBases * percentile);
        long cumulativeTotal = 0;

        for(int i = 0; i < baseQualFrequency.length; ++i)
        {
            cumulativeTotal += baseQualFrequency[i];

            if(cumulativeTotal >= percentileEntry)
                return (byte) i;
        }

        return LOW_BASE_QUAL_THRESHOLD;
    }

    public static Map<HlaGene, int[]> calculateGeneCoverage(final Iterable<Fragment> fragments)
    {
        final Map<HlaGene, int[]> geneBaseDepth = Maps.newHashMap();

        for(HlaGene geneName : GENE_CACHE.GeneNames)
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
