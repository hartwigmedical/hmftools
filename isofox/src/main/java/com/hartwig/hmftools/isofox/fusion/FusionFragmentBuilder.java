package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.hasRealignableSoftClip;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isRealignedFragmentCandidate;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FusionFragmentBuilder
{
    /*
    Primary fusion fragments:
    - 3 reads including supplementary
        - may be split across gene collections, or a local INV, DUP or a very long DEL
        - covers all SV types - DEL, DUP, BND and INV
    - 2 reads in same gene collection with DUP orientations - should still have a supplementary
    - 2 reads in same gene collection with N-split and clearly involving 2 genes (eg at Known/Canonical sites), overlapped by longer intronic gene
        - very hard to know at time of BAM processing since could always be intronic..

    Secondary fusion fragments:
    - 2 reads forming a local DEL without supplementary where 1 or more reads ends outside the current gene collection
        - has an N-split for the actual DEL
        - has a pair of soft-clipped facing reads
    - 2 discordant reads if not supporting a fusion and in a known pair

    Supporting fragments:
    - 2 reads, one in each gene collection - discordant reads supporting a fusion
    - 2 reads in same gene collection with soft-clipping, potential realignable bases

     */

    public static void setFragmentProperties(final FusionFragment fragment)
    {
        /* Set the following propeties for fragments:
        - chromosomes
        - junction properties if clear - position, orientation
        - gene collection
        - fragment type - junction, discordant or unknown (realigned will be set later)

         Scenarios for fragment junctions and read positions
        - interchromosomal - typically with a supplementary but otherwise dual soft-clipping (possibly discordant?)
        - single-gene inversion with supplementary alignment
        - inter-gene INV, DUP or DEL with supplementary alignment
        - split-read inter-gene either within the same gene collection or spanning them

        Non-2-junction fragments:
        - Single-junction - realignable single junction fragments - may or may not be single-GC
        - Discordant - currently must span gene collections to be considered
        */

        // set chromosomes - either 1 or 2, with the lower set into the start position
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);

        for(final FusionRead read : fragment.reads())
        {
            if(!chromosomes.contains(read.Chromosome))
                chromosomes.add(read.Chromosome);
        }

        if(chromosomes.size() == 1)
        {
            fragment.chromosomes()[SE_START] = fragment.chromosomes()[SE_END] = chromosomes.get(0);
        }
        else
        {
            if(lowerChromosome(chromosomes.get(0), chromosomes.get(1)))
            {
                fragment.chromosomes()[SE_START] = chromosomes.get(0);
                fragment.chromosomes()[SE_END] = chromosomes.get(1);
            }
            else
            {
                fragment.chromosomes()[SE_START] = chromosomes.get(1);
                fragment.chromosomes()[SE_END] = chromosomes.get(0);
            }
        }

        // first handle reads with supplementary alignment since these should be most clear
        if(fragment.hasSuppAlignment())
        {
            setSuppAlignJunctionData(fragment);
            return;
        }

        // identify the properties of the junction and use it to set other properties of the fragment
        if(fragment.isSingleChromosome())
        {
            final FusionRead splitRead = fragment.reads().stream()
                    .filter(x -> x.ContainsSplit)
                    .filter(x -> x.spansGeneCollections() || x.HasInterGeneSplit)
                    .findFirst().orElse(null);

            if(splitRead != null)
            {
                // set the junction data around the spanning N-split
                setSplitReadJunctionData(fragment, splitRead);
                return;
            }

            // set single junction info for candidate realignable fragments
            if(fragment.reads().size() == 2
            && fragment.reads().stream().anyMatch(x -> isRealignedFragmentCandidate(x))
            && fragment.reads().stream().noneMatch(x -> x.spansGeneCollections()))
            {
                FusionRead read1 = fragment.reads().get(0);
                FusionRead read2 = fragment.reads().get(1);
                if(read1.GeneCollections[0] == read2.GeneCollections[1])
                {
                    setSingleSoftClipJunctionData(fragment);
                    return;
                }
            }
        }

        // set gene collection and orientation for discordant reads
        setNonJunctionData(fragment);
    }

    private static void setSuppAlignJunctionData(final FusionFragment fragment)
    {
        int[] junctionPositions = { -1, -1 };
        byte[] junctionOrientations = { 0, 0 };
        int[] geneCollections = { -1, -1 };
        String[] chromosomes = {"", ""};

        int posIndex = 0;

        for(FusionRead read : fragment.reads())
        {
            if(!read.HasSuppAlignment)
                continue;

            chromosomes[posIndex] = read.Chromosome;

            if(read.isLongestSoftClip(SE_START))
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_START);
                junctionOrientations[posIndex] = ORIENT_REV;
                geneCollections[posIndex] = read.GeneCollections[SE_START];
            }
            else
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_END);
                junctionOrientations[posIndex] = ORIENT_FWD;
                geneCollections[posIndex] = read.GeneCollections[SE_END];
            }

            ++posIndex;

            if(posIndex >= 2)
                break;
        }

        // set orientations and lower/upper index based on position
        if(junctionPositions[0] < 0 || junctionPositions[1] < 0)
            return;

        int lowerIndex;

        if(chromosomes[SE_START].equals(chromosomes[SE_END]))
            lowerIndex = junctionPositions[0] <= junctionPositions[1] ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes[SE_START], chromosomes[SE_END]) ? 0 : 1;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            fragment.junctionPositions()[se] = junctionPositions[index];
            fragment.junctionOrientations()[se] = junctionOrientations[index];
            fragment.orientations()[se] = junctionOrientations[index];
            fragment.geneCollections()[se] = geneCollections[index];
        }

        fragment.setType(MATCHED_JUNCTION);
    }

    private static void setSplitReadJunctionData(final FusionFragment fragment, final FusionRead splitRead)
    {
        // set the junction data around the spanning N-split
        final int[] splitJunction = findSplitReadJunction(splitRead);

        if(splitJunction != null)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                fragment.geneCollections()[se] = splitRead.GeneCollections[se];
                fragment.junctionPositions()[se] = splitJunction[se];
            }

            fragment.junctionOrientations()[SE_START] = fragment.orientations()[SE_START] = ORIENT_FWD;
            fragment.junctionOrientations()[SE_END] = fragment.orientations()[SE_END] = ORIENT_REV;
            fragment.setType(MATCHED_JUNCTION);
        }
    }

    private static void setSingleSoftClipJunctionData(final FusionFragment fragment)
    {
        // 2 gene collections are involved or a candidate realignable fragment in one collection
        FusionRead realignRead = null;
        int maxScLength = 0;
        int scSide = 0;

        for(FusionRead read : fragment.reads())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(!hasRealignableSoftClip(read, se, true))
                    continue;

                if(read.SoftClipLengths[se] > maxScLength)
                {
                    maxScLength = read.SoftClipLengths[se];
                    scSide = se;
                    realignRead = read;
                }
            }
        }

        if(realignRead == null)
            return;

        // by convention store the junction in the start slot since the other potential junction is unknown at this stage
        fragment.junctionPositions()[SE_START] = realignRead.getCoordsBoundary(scSide);
        fragment.junctionOrientations()[SE_START] = scSide == SE_START ? ORIENT_REV : ORIENT_FWD;
        fragment.orientations()[SE_START] = fragment.junctionOrientations()[SE_START];

        fragment.geneCollections()[SE_START] = fragment.geneCollections()[SE_END] = realignRead.GeneCollections[scSide];
        fragment.setType(REALIGN_CANDIDATE);
    }

    private static void setNonJunctionData(final FusionFragment fragment)
    {
        // set gene collections from the reads without any knowledge of which junctions they may support
        final List<String> chrGeneCollections = Lists.newArrayListWithCapacity(2);
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);
        final Map<String,Integer> positions = Maps.newHashMapWithExpectedSize(2);
        final List<Integer> geneCollections = Lists.newArrayListWithCapacity(2);
        final Map<String,FusionRead> reads = Maps.newHashMapWithExpectedSize(2);

        final Map<String,List<FusionRead>> readGroups = Maps.newHashMapWithExpectedSize(2);

        for(final FusionRead read : fragment.reads())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                // skip the upper side since it is likely to be the non-genic region of the next gene collection
                if(!read.spansGeneCollections() && se == SE_END)
                    continue;

                final String chrGeneId = formLocation(read.Chromosome, read.GeneCollections[se], true); // genic status ignored for group determination

                List<FusionRead> readGroup = readGroups.get(chrGeneId);

                if(readGroup == null)
                {
                    readGroups.put(chrGeneId, Lists.newArrayList(read));

                    chrGeneCollections.add(chrGeneId);
                    chromosomes.add(read.Chromosome);
                    geneCollections.add(read.GeneCollections[se]);
                    positions.put(chrGeneId, read.getCoordsBoundary(se)); // no overlap in gene collections so doesn't matter which position is used
                    reads.put(chrGeneId, read);
                }
                else
                {
                    if(read.getCoordsBoundary(se) < positions.get(chrGeneId))
                    {
                        positions.put(chrGeneId, min(positions.get(chrGeneId), read.getCoordsBoundary(se)));
                        reads.put(chrGeneId, read);
                    }

                    readGroup.add(read);
                }
            }
        }

        // now set start and end chromosome, orientation and gene collections from amongst the reads in each group
        if(readGroups.size() == 1)
        {
            fragment.geneCollections()[SE_START] = fragment.geneCollections()[SE_END] = geneCollections.get(0);

            // orientation could be set based on the orientations and positions of the reads.. do this when the junction data is set
            fragment.orientations()[SE_START] = ORIENT_FWD;
            fragment.orientations()[SE_END] = ORIENT_REV;
            return;
        }

        // find the lower chromosome and position
        int lowerIndex;

        if(chromosomes.get(0).equals(chromosomes.get(1)))
            lowerIndex = positions.get(chrGeneCollections.get(0)) <= positions.get(chrGeneCollections.get(1)) ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes.get(0), chromosomes.get(1)) ? 0 : 1;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);

            fragment.geneCollections()[se] = geneCollections.get(index);
            fragment.orientations()[se] = reads.get(chrGeneId).Orientation;
        }

        fragment.setType(DISCORDANT);

        if(readGroups.size() != 2 || fragment.reads().size() != 2 || !fragment.isSingleChromosome())
            return;

        // check for evidence of a fusion junction from soft-clipped reads but without supplementary alignment data
        int[] scPositions = {-1, -1};
        byte[] scOrientations = {0, 0};

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);
            FusionRead read = reads.get(chrGeneId);

            int requiredScSide = switchIndex(se);

            if(!read.isSoftClipped(requiredScSide))
                break;

            if(requiredScSide == SE_START)
            {
                scPositions[se] = read.getCoordsBoundary(SE_START);
                scOrientations[se] = ORIENT_REV;
            }
            else
            {
                scPositions[se] = read.getCoordsBoundary(SE_END);
                scOrientations[se] = ORIENT_FWD;
            }
        }

        if(scPositions[SE_START] > 0 && scPositions[SE_START] < scPositions[SE_END])
        {
            // this will now be a candidate for a local fusion from a discordant read pair
            fragment.setType(DISCORDANT_JUNCTION);

            for(int se = SE_START; se <= SE_END; ++se)
            {
                int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
                fragment.junctionPositions()[se] = scPositions[index];
                fragment.junctionOrientations()[se] = scOrientations[index];
            }
        }
    }
}
