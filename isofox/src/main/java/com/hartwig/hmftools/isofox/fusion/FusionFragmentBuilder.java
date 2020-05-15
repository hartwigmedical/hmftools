package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.NEG_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.hasRealignableSoftClip;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isInversion;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isRealignedFragmentCandidate;
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

        // dismiss any read with more than 2

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
        return reads.stream().anyMatch(x -> x.hasSuppAlignment());
    }

    public static boolean spansGeneCollections(final List<ReadRecord> reads)
    {
        if(reads.stream().anyMatch(x -> x.spansGeneCollections()))
            return true;

        int gc = reads.get(0).getGeneCollectons()[SE_START];
        final String chr = reads.get(0).Chromosome;

        for(int i = 1; i < reads.size(); ++i)
        {
            if(!chr.equals(reads.get(i).Chromosome) || reads.get(i).getGeneCollectons()[SE_START] != gc)
                return true;
        }

        return false;
    }

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

        for(final ReadRecord read : fragment.reads())
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
        if(fragment.isSingleChromosome()
        && fragment.reads().stream().anyMatch(x -> x.containsSplit() && (x.spansGeneCollections() || x.hasInterGeneSplit())))
        {
            // set the junction data around the spanning N-split
            setSplitReadJunctionData(fragment);
            return;
        }

        // set single junction info for candidate realignable reads
        if(fragment.reads().stream().anyMatch(x -> isRealignedFragmentCandidate(x)))
        {
            setSingleSoftClipJunctionData(fragment);
            return;
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

        for(ReadRecord read : fragment.reads())
        {
            if (!read.hasSuppAlignment())
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

            chromosomes[posIndex] = read.Chromosome;

            if (useLeft)
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_START);
                junctionOrientations[posIndex] = NEG_ORIENT;
                geneCollections[posIndex] = read.getGeneCollectons()[SE_START];
            }
            else
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_END);
                junctionOrientations[posIndex] = POS_ORIENT;
                geneCollections[posIndex] = read.getGeneCollectons()[SE_END];
            }

            ++posIndex;

            if(posIndex >= 2)
                break;
        }

        // set orientations and lower/upper index based on position
        if(junctionPositions[0] < 0 || junctionPositions[1] < 0)
            return;

        int lowerIndex;

        if (chromosomes[SE_START].equals(chromosomes[SE_END]))
            lowerIndex = junctionPositions[0] <= junctionPositions[1] ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes[SE_START], chromosomes[SE_END]) ? 0 : 1;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            fragment.junctionPositions()[se] = junctionPositions[index];
            fragment.junctionOrientations()[se] = junctionOrientations[index];
            fragment.orientations()[se] = junctionOrientations[index];
            fragment.geneCollections()[se] = geneCollections[index];
        }

        fragment.setType(MATCHED_JUNCTION);
    }

    private static void setSplitReadJunctionData(final FusionFragment fragment)
    {
        // set the junction data around the spanning N-split
        final ReadRecord splitRead = fragment.reads().stream()
                .filter(x -> x.spansGeneCollections() || x.hasInterGeneSplit())
                .findFirst().orElse(null);

        final int[] splitJunction = findSplitReadJunction(splitRead);

        if(splitJunction != null)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                fragment.geneCollections()[se] = splitRead.getGeneCollectons()[se];
                fragment.junctionPositions()[se] = splitJunction[se];
            }

            fragment.junctionOrientations()[SE_START] = fragment.orientations()[SE_START] = POS_ORIENT;
            fragment.junctionOrientations()[SE_END] = fragment.orientations()[SE_END] = NEG_ORIENT;
            fragment.setType(MATCHED_JUNCTION);
        }
    }

    private static void setSingleSoftClipJunctionData(final FusionFragment fragment)
    {
        // 2 gene collections are involved or a candidate realignable fragment in one collection
        ReadRecord realignRead = null;
        int maxScLength = 0;
        int scSide = 0;

        for(ReadRecord read : fragment.reads())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(!hasRealignableSoftClip(read, se, true))
                    continue;

                int scLength = se == SE_START ? read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

                if(scLength > maxScLength)
                {
                    maxScLength = scLength;
                    scSide = se;
                    realignRead = read;
                }
            }
        }

        if(realignRead == null)
            return;

        // by convention store the junction in the start slot since the other potential junction is unknown at this stage
        fragment.junctionPositions()[SE_START] = realignRead.getCoordsBoundary(scSide);
        fragment.junctionOrientations()[SE_START] = scSide == SE_START ? NEG_ORIENT : POS_ORIENT;
        fragment.orientations()[SE_START] = fragment.junctionOrientations()[SE_START];

        fragment.geneCollections()[SE_START] = fragment.geneCollections()[SE_END] = realignRead.getGeneCollectons()[scSide];
        fragment.setType(REALIGN_CANDIDATE);
    }

    private static void setNonJunctionData(final FusionFragment fragment)
    {
        // set gene collections from the reads without any knowledge of which junctions they may support
        final List<String> chrGeneCollections = Lists.newArrayListWithCapacity(2);
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);
        final Map<String,Integer> positions = Maps.newHashMapWithExpectedSize(2);
        final List<Integer> geneCollections = Lists.newArrayListWithCapacity(2);
        final Map<String,ReadRecord> reads = Maps.newHashMapWithExpectedSize(2);

        final Map<String,List<ReadRecord>> readGroups = Maps.newHashMapWithExpectedSize(2);

        for(final ReadRecord read : fragment.reads())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                // skip the upper side since it is likely to be the non-genic region of the next gene collection
                if(!read.spansGeneCollections() && se == SE_END)
                    continue;

                final String chrGeneId = formLocation(read.Chromosome, read.getGeneCollectons()[se], true); // genic status ignored for group determination

                List<ReadRecord> readGroup = readGroups.get(chrGeneId);

                if (readGroup == null)
                {
                    readGroups.put(chrGeneId, Lists.newArrayList(read));

                    chrGeneCollections.add(chrGeneId);
                    chromosomes.add(read.Chromosome);
                    geneCollections.add(read.getGeneCollectons()[se]);
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
            fragment.orientations()[SE_START] = POS_ORIENT;
            fragment.orientations()[SE_END] = NEG_ORIENT;

            return;
        }

        // find the lower chromosome and position
        int lowerIndex;

        if (chromosomes.get(0).equals(chromosomes.get(1)))
            lowerIndex = positions.get(chrGeneCollections.get(0)) <= positions.get(chrGeneCollections.get(1)) ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes.get(0), chromosomes.get(1)) ? 0 : 1;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);

            fragment.geneCollections()[se] = geneCollections.get(index);
            fragment.orientations()[se] = reads.get(chrGeneId).orientation();
        }

        fragment.setType(DISCORDANT);

        if(readGroups.size() != 2 || fragment.reads().size() != 2 || !fragment.isSingleChromosome())
            return;

        // check for evidence of a fusion junction from soft-clipped reads but without supplementary alignment data
        int[] scPositions = {-1, -1};
        byte[] scOrientations = {0, 0};

        for (int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);
            ReadRecord read = reads.get(chrGeneId);

            int requiredScSide = switchIndex(se);

            if(!read.isSoftClipped(requiredScSide))
                break;

            int scLength = requiredScSide == SE_START ?
                    read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

            if(scLength < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH)
                break;

            if(requiredScSide == SE_START)
            {
                scPositions[se] = read.getCoordsBoundary(SE_START);
                scOrientations[se] = -1;
            }
            else
            {
                scPositions[se] = read.getCoordsBoundary(SE_END);
                scOrientations[se] = 1;
            }
        }

        if(scPositions[SE_START] > 0 && scPositions[SE_END] > 0)
        {
            // this will now be a candidate for a local fusion
            fragment.setType(MATCHED_JUNCTION);

            for (int se = SE_START; se <= SE_END; ++se)
            {
                int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
                fragment.junctionPositions()[se] = scPositions[index];
                fragment.junctionOrientations()[se] = scOrientations[index];
            }
        }
    }


    //////////////
    // methods to be deprecated
    /*
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
                        .filter(x -> !x.hasSuppAlignment()).findFirst().map(x -> x.orientation()).orElse((byte)0);
            else
                fragment.orientations()[se] = readGroup.get(0).orientation();
        }
    }

    public static void setJunctionData(final FusionFragment fragment)
    {
        if(fragment.reads().stream().anyMatch(x -> x.spansGeneCollections()))
        {
            // set the junction data around the spanning N-split
            setSplitReadJunctionData(fragment);
            return;
        }

        if(fragment.hasSuppAlignment() && fragment.isSingleGene())
        {
            if(!isInversion(fragment.reads()))
            {
                ISF_LOGGER.warn("read({}) has supp in single gene but not INV", fragment.reads().get(0));
            }

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
                    if(rgRead.hasSuppAlignment())
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

    public static String[] setLocationIds(final FusionFragment fragment)
    {
        final String[] locationIds = new String[SE_PAIR];

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(fragment.isSingleGene() && se == SE_END)
                continue;

            final List<ReadRecord> readsOnSide = fragment.readsByLocation(se);

            if (fragment.junctionPositions()[se] > 0)
            {
                final int seIndex = se;
                int closestDistance = -1;
                int juncGeneCollection = 0;

                for(ReadRecord read : readsOnSide)
                {
                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        int distanceFromJunc = abs(read.getCoordsBoundary(se2) - fragment.junctionPositions()[se]);
                        if(closestDistance == -1 || distanceFromJunc < closestDistance)
                        {
                            closestDistance = distanceFromJunc;
                            juncGeneCollection = read.getGeneCollectons()[se2];
                        }
                    }
                }

                locationIds[se] = formLocation(fragment.chromosomes()[se], juncGeneCollection, true);
            }
            else
            {
                // for the discordant case, where a read spans gene collectons, take the lower read end since it will be genic
                if(!readsOnSide.isEmpty())
                {
                    ReadRecord read = readsOnSide.get(0);
                    locationIds[se] = formLocation(read.Chromosome, read.getGeneCollectons()[SE_START], true);
                }
            }
        }

        return locationIds;
    }

    private static void setSingleGeneSuppAlignJunctionData(final FusionFragment fragment)
    {
        int[] junctionPositions = { -1, -1 };
        byte[] junctionOrientations = { 0, 0 };
        int[] geneCollections = { -1, -1 };

        int posIndex = 0;

        for(ReadRecord read : fragment.reads())
        {
            if (!read.hasSuppAlignment())
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
                geneCollections[posIndex] = read.getGeneCollectons()[SE_START];
            }
            else
            {
                junctionPositions[posIndex] = read.getCoordsBoundary(SE_END);
                junctionOrientations[posIndex] = 1;
                geneCollections[posIndex] = read.getGeneCollectons()[SE_END];
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

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? startIndex : endIndex;
            fragment.junctionPositions()[se] = junctionPositions[index];
            fragment.junctionOrientations()[se] = junctionOrientations[index];
            fragment.orientations()[se] = junctionOrientations[index];
            fragment.geneCollections()[se] = geneCollections[index];
        }

        fragment.setType(MATCHED_JUNCTION);
    }

     */

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
