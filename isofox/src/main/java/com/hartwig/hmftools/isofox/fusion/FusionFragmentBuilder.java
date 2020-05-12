package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.lowerChromosome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.common.ReadRecord;

import htsjdk.samtools.CigarOperator;

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
        - has soft-clipped facing reads
    - 2 discordant reads if not supporting a fusion and in a known pair

    Supporting fragments:
    - 2 reads, one in each gene collection - discordant reads supporting a fusion
    - 2 reads in same gene collection with soft-clipping, potential realignable bases

     */

    public static boolean isValidFragment(final List<ReadRecord> reads)
    {
        if(reads.size() <= 1)
            return false;

        if(hasSuppAlignment(reads))
        {
            return (reads.size() == 3);
        }

        // no further conditions

        /*
        Set<String> chrGeneSet = Sets.newHashSetWithExpectedSize(3);

        for(ReadRecord read : reads)
        {
            chrGeneSet.add(formLocation(read.Chromosome, read.getGeneCollecton()));
            if(chrGeneSet.size() == 3)
                return false;
        }

        if(chrGeneSet.size() == 2)
            return true;

        if(reads.stream().anyMatch(x -> isRealignedFragmentCandidate(x)))
            return true;

        return false;

        */

        return (reads.size() == 2);
    }

    public static boolean hasSuppAlignment(final List<ReadRecord> reads)
    {
        return reads.stream().anyMatch(x -> x.getSuppAlignment() != null);
    }

    public static void setLocationData(final FusionFragment fragment)
    {
        // divide the reads into groups by chromosome and gene collection - so reads may span both
        final List<String> chrGeneCollections = Lists.newArrayListWithCapacity(2);
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);
        final List<Integer> positions = Lists.newArrayListWithCapacity(2);
        final List<Integer> geneCollections = Lists.newArrayListWithCapacity(2);
        final List<Boolean> inGenicRegions = Lists.newArrayListWithCapacity(2);

        final Map<String,List<ReadRecord>> readGroups = Maps.newHashMapWithExpectedSize(2);

        for(final ReadRecord read : fragment.reads())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(!read.spansGeneCollections() && se == SE_END)
                    continue;

                int geneId = read.getGeneCollectons()[se];

                final String chrGeneId = formLocation(read.Chromosome, geneId, read.getIsGenicRegion()[se]);

                List<ReadRecord> readGroup = readGroups.get(chrGeneId);

                if (readGroup == null)
                {
                    readGroups.put(chrGeneId, Lists.newArrayList(read));

                    chrGeneCollections.add(chrGeneId);
                    chromosomes.add(read.Chromosome);
                    geneCollections.add(read.getGeneCollectons()[se]);
                    inGenicRegions.add(read.getIsGenicRegion()[se]);
                    positions.add(read.getCoordsBoundary(se)); // no overlap in gene collections so doesn't matter which position is used
                }
                else
                {
                    readGroup.add(read);
                }
            }
        }

        // now set start and end chromosome, orientation and gene collections from amongst the reads in each group
        if(readGroups.size() == 1)
        {
            fragment.chromosomes()[SE_START] = fragment.chromosomes()[SE_END] = chromosomes.get(0);
            fragment.geneCollections()[SE_START] = fragment.geneCollections()[SE_END] = geneCollections.get(0);
            fragment.inGenicRegions()[SE_START] = fragment.inGenicRegions()[SE_END] = inGenicRegions.get(0);
            fragment.readsByLocation(SE_START).addAll(fragment.reads());

            // orientation could be set based on the orientations and positions of the reads.. do this when the junction data is set
            fragment.orientations()[SE_START] = fragment.orientations()[SE_END] = fragment.reads().get(0).orientation();

            return;
        }

        // find the lower chromosome and position
        int lowerIndex;

        if (chromosomes.get(0).equals(chromosomes.get(1)))
            lowerIndex = positions.get(0) <= positions.get(1) ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes.get(0), chromosomes.get(1)) ? 0 : 1;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);

            final List<ReadRecord> readGroup = readGroups.get(chrGeneId);

            fragment.readsByLocation(se).addAll(readGroup);
            fragment.chromosomes()[se] = chromosomes.get(index);
            fragment.geneCollections()[se] = geneCollections.get(index);
            fragment.inGenicRegions()[se] = inGenicRegions.get(index);

            if(readGroup.size() == 2)
                fragment.orientations()[se] = readGroup.stream()
                        .filter(x -> x.getSuppAlignment() == null).findFirst().map(x -> x.orientation()).orElse((byte)0);
            else
                fragment.orientations()[se] = readGroup.get(0).orientation();
        }
    }

    private static void setSplitReadJunctionData(final FusionFragment fragment)
    {
        // set the junction data around the spanning N-split
        final ReadRecord splitRead = fragment.reads().stream().filter(x -> x.spansGeneCollections()).findFirst().orElse(null);

        final int[] splitJunction = findSplitReadJunction(splitRead);

        if(splitJunction != null)
        {
            fragment.junctionOrientations()[SE_START] = 1;
            fragment.junctionOrientations()[SE_END] = -1;
            fragment.junctionPositions()[SE_START] = splitJunction[SE_START];
            fragment.junctionPositions()[SE_END] = splitJunction[SE_END];
            fragment.setType(MATCHED_JUNCTION);
            return;
        }
    }

    private static void setSingleGeneSuppAlignJunctionData(final FusionFragment fragment)
    {
        int[] junctionPositions = { -1, -1 };
        byte[] junctionOrientations = { 0, 0 };

        int posIndex = 0;

        for(ReadRecord read : fragment.readsByLocation(SE_START))
        {
            if (read.getSuppAlignment() == null)
                continue;

            int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
            int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

            boolean useLeft = false;

            if (scLeft > 0 && scRight > 0)
            {
                if (scLeft >= scRight)
                {
                    useLeft = true;
                }
                else if (scRight > scLeft)
                {
                    useLeft = false;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                useLeft = scLeft > 0;
            }

            if (useLeft)
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_START);
                junctionOrientations[posIndex] = -1;
            }
            else
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_END);
                junctionOrientations[posIndex] = 1;
            }

            ++posIndex;

            if(posIndex >= 2)
                break;
        }

        // set orientations and lower/upper index based on position
        if(junctionPositions[0] < 0 || junctionPositions[1] < 0)
            return;

        int startIndex = junctionPositions[0] < junctionPositions[1] ? 0 : 1;
        int endIndex = switchIndex(startIndex);

        fragment.junctionPositions()[SE_START] = junctionPositions[startIndex];
        fragment.junctionPositions()[SE_END] = junctionPositions[endIndex];
        fragment.junctionOrientations()[SE_START] = junctionOrientations[startIndex];
        fragment.junctionOrientations()[SE_END] = junctionOrientations[endIndex];
        fragment.orientations()[SE_START] = junctionOrientations[startIndex];
        fragment.orientations()[SE_END] = junctionOrientations[endIndex];

        fragment.setType(MATCHED_JUNCTION);
    }

    public static void setJunctionData(final FusionFragment fragment)
    {
        /* junctions can be identified by a few methods:
            - supplementary alignment info - within the same gene collection (eg a DUP, INV or very long DEL) or not (anything)
            - 2 reads forming a DEL with facing, matching soft-clipped bases in different gene collections
            - a N-split DEL within the same read linking 2 genes (currently must be in 2 diff gene collections)
            - a candidate realignment fragment where only 1 junction is supported
        */

        if(fragment.reads().stream().anyMatch(x -> x.spansGeneCollections()))
        {
            // set the junction data around the spanning N-split
            setSplitReadJunctionData(fragment);
            return;
        }

        if(fragment.hasSuppAlignment() && fragment.isSingleGene())
        {
            setSingleGeneSuppAlignJunctionData(fragment);
            return;
        }

        // 2 gene collections are involved or a candidate realignable fragment in one collection

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(fragment.isSingleGene() && se == SE_END)
                continue;

            final List<ReadRecord> readGroup = fragment.readsByLocation(se);

            if(readGroup.isEmpty())
                continue;

            ReadRecord read = null;
            int maxScLength = 0;

            for(ReadRecord rgRead : readGroup)
            {
                if(!rgRead.Cigar.containsOperator(CigarOperator.S))
                    continue;

                if(fragment.hasSuppAlignment())
                {
                    if(rgRead.getSuppAlignment() != null)
                        read = rgRead;
                    else
                        continue;
                }
                else
                {
                    int scLeft = rgRead.isSoftClipped(SE_START) ? rgRead.Cigar.getFirstCigarElement().getLength() : 0;
                    int scRight = rgRead.isSoftClipped(SE_END) ? rgRead.Cigar.getLastCigarElement().getLength() : 0;

                    if(max(scLeft, scRight) > maxScLength)
                    {
                        maxScLength = max(scLeft, scRight);
                        read = rgRead;
                    }
                }
            }

            if(read == null)
                continue;

            int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
            int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

            boolean useLeft = false;

            if(scLeft > 0 && scRight > 0)
            {
                // should be very unlikely since implies a very short exon and even then would expect it to be mapped
                if(scLeft >= scRight)
                {
                    useLeft = true;
                }
                else if(scRight > scLeft)
                {
                    useLeft = false;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                useLeft = scLeft > 0;
            }

            if(useLeft)
            {
                fragment.junctionPositions()[se] = read.getCoordsBoundary(SE_START);
                fragment.junctionOrientations()[se] = -1;
            }
            else
            {
                fragment.junctionPositions()[se] = read.getCoordsBoundary(SE_END);
                fragment.junctionOrientations()[se] = 1;
            }
        }

        if(fragment.junctionPositions()[SE_START] > 0 && fragment.junctionPositions()[SE_END] > 0)
        {
            fragment.setType(MATCHED_JUNCTION);
        }
        else if(!fragment.isSingleGene())
        {
            fragment.setType(DISCORDANT);
        }
    }

    /*  Junction bases will now be set from the ref genome, not the read, in case they contain SNVs and cause invalid realignments

        int baseLength = min(10, read.Length);

        if(mJunctionOrientations[se] == 1)
        {
            int scLength = read.Cigar.getLastCigarElement().getLength();
            int readEndPos = read.Length - scLength;
            mJunctionBases[se] = read.ReadBases.substring(readEndPos - baseLength, readEndPos);
            // softClipBases[se] = read.ReadBases.substring(readEndPos, readEndPos + baseLength);
        }
        else
        {
            int scLength = read.Cigar.getFirstCigarElement().getLength();
            int readStartPos = scLength;
            mJunctionBases[se] = read.ReadBases.substring(readStartPos, readStartPos + baseLength);
            // softClipBases[se] = read.ReadBases.substring(readStartPos - baseLength, readStartPos);
        }

    if(!mJunctionBases[SE_START].isEmpty() && !mJunctionBases[SE_END].isEmpty())
    {
        if(mJunctionBases[SE_START].equals(softClipBases[SE_END]) && mJunctionBases[SE_END].equals(softClipBases[SE_START]))
        {
            mJunctionBasesMatched = true;
        }
    }
    */


}
