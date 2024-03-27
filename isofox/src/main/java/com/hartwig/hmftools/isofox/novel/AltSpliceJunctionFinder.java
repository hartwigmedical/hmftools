package com.hartwig.hmftools.isofox.novel;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.EXONIC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.SPLICE_JUNC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.CIRCULAR;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.EXON_INTRON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.INTRONIC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.MIXED_TRANS;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_INTRON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.SKIPPED_EXONS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransMatchType;

import htsjdk.samtools.CigarOperator;

public class AltSpliceJunctionFinder
{
    private final boolean mEnabled;
    private final IsofoxConfig mConfig;

    private final List<AltSpliceJunction> mAltSpliceJunctions;
    private final BufferedWriter mWriter;

    private GeneCollection mGenes;

    public AltSpliceJunctionFinder(final IsofoxConfig config, final BufferedWriter writer)
    {
        mEnabled = config.runFunction(ALT_SPLICE_JUNCTIONS);
        mConfig = config;
        mAltSpliceJunctions = Lists.newArrayList();
        mWriter = writer;
        mGenes = null;
    }

    public boolean enabled() { return mEnabled; }
    public List<AltSpliceJunction> getAltSpliceJunctions() { return mAltSpliceJunctions; }

    public void setGeneData(final GeneCollection genes)
    {
        mGenes = genes;
        mAltSpliceJunctions.clear();
    }

    public void evaluateFragmentReads(
            final List<GeneReadData> genes, final ReadRecord read1, final ReadRecord read2, final List<Integer> relatedTransIds)
    {
        if(!mEnabled)
            return;

        if(read1.isDuplicate() || read2.isDuplicate() || genes.isEmpty() || read1.isMultiMapped() || read2.isMultiMapped())
            return;

        // exclude SJs too far outside known transcripts
        int[] geneBounds = new int[] {
                genes.stream().mapToInt(x -> x.GeneData.GeneStart).min().orElse(0) - MAX_NOVEL_SJ_DISTANCE,
                genes.stream().mapToInt(x -> x.GeneData.GeneStart).max().orElse(0) + MAX_NOVEL_SJ_DISTANCE };

        if(!positionsWithin(read1.PosStart, read1.PosEnd, geneBounds[SE_START], geneBounds[SE_END])
        || !positionsWithin(read2.PosStart, read2.PosEnd, geneBounds[SE_START], geneBounds[SE_END]))
        {
            return;
        }

        // at least one of the reads must fall within a gene
        final List<GeneReadData> candidateGenes = genes.stream()
                .filter(x -> positionsWithin(read1.PosStart, read1.PosEnd, x.GeneData.GeneStart,x.GeneData.GeneEnd)
                        || positionsWithin(read2.PosStart, read2.PosEnd, x.GeneData.GeneStart,x.GeneData.GeneEnd))
                .collect(Collectors.toList());

        if(candidateGenes.isEmpty())
            return;

        AltSpliceJunction firstAltSJ = null;

        if(isCandidateCircular(read1, read2))
        {
            firstAltSJ = registerAltSpliceJunction(candidateGenes, read1, read2, relatedTransIds);

            if(firstAltSJ == null)
                return;
        }
        else
        {
            if(AltSpliceJunctionFinder.isCandidate(read1))
            {
                firstAltSJ = registerAltSpliceJunction(candidateGenes, read1, relatedTransIds);
            }

            AltSpliceJunction secondAltSJ = null;

            if(AltSpliceJunctionFinder.isCandidate(read2))
            {
                // avoid double-counting overlapping reads
                if(firstAltSJ != null && positionsOverlap(read1.PosStart, read1.PosEnd, read2.PosStart, read2.PosEnd))
                    return;

                secondAltSJ = registerAltSpliceJunction(candidateGenes, read2, relatedTransIds);
            }

            if(firstAltSJ != null && secondAltSJ != null)
            {
                checkNovelExon(firstAltSJ, secondAltSJ);
            }
        }
    }

    private static boolean isCandidate(final ReadRecord read)
    {
        if(!read.Cigar.containsOperator(CigarOperator.N))
            return false;

        if(read.getTranscriptClassifications().values().contains(TransMatchType.SPLICE_JUNCTION))
            return false;

        if(read.getMappedRegionCoords().size() == 1)
            return false;

        return true;
    }

    private static boolean isCandidateCircular(final ReadRecord read1, final ReadRecord read2)
    {
        return read1.hasSuppAlignment() && read2.hasSuppAlignment();
    }

    public boolean junctionMatchesKnownSpliceJunction(final List<GeneReadData> genes, final int[] spliceJunction, final List<Integer> transIds)
    {
        for(GeneReadData gene : genes)
        {
            for (TranscriptData transData : gene.getTranscripts())
            {
                if (!transIds.contains(transData.TransId))
                    continue;

                for (int i = 0; i < transData.exons().size() - 1; ++i)
                {
                    if (transData.exons().get(i).End == spliceJunction[SE_START]
                    && transData.exons().get(i + 1).Start == spliceJunction[SE_END])
                        return true;
                }
            }
        }

        return false;
    }

    public AltSpliceJunction createFromRead(final ReadRecord read, final List<Integer> relatedTransIds)
    {
        // related transcripts will any of those where either read covers one or more of its exons
        int[] spliceJunction = new int[SE_PAIR];

        // find the novel splice junction, and all associated transcripts
        final List<int[]> mappedCoords = read.getMappedRegionCoords();

        if(read.inferredCoordAdded(true))
        {
            spliceJunction[SE_START] = mappedCoords.get(1)[SE_END];
            spliceJunction[SE_END] = mappedCoords.get(2)[SE_START];
        }
        else
        {
            spliceJunction[SE_START] = mappedCoords.get(0)[SE_END];
            spliceJunction[SE_END] = mappedCoords.get(1)[SE_START];
        }

        checkJunctionHomology(spliceJunction, read.getMappedRegions());

        List<RegionReadData> sjStartRegions = Lists.newArrayList(); // transcript regions with an exon matching the start of the alt SJ
        List<RegionReadData> sjEndRegions = Lists.newArrayList();
        final AltSpliceJunctionContext[] regionContexts = { AltSpliceJunctionContext.UNKNOWN, AltSpliceJunctionContext.UNKNOWN };

        classifyRegions(read, spliceJunction, sjStartRegions, sjEndRegions, regionContexts);

        AltSpliceJunctionType sjType = classifySpliceJunction(relatedTransIds, sjStartRegions, sjEndRegions, regionContexts);

        AltSpliceJunction altSplicJunction = new AltSpliceJunction(
                mGenes.chromosome(), spliceJunction, sjType, read.Id, regionContexts, sjStartRegions, sjEndRegions);

        altSplicJunction.setCandidateTranscripts(read.getMappedRegions().keySet().stream().collect(Collectors.toList()));

        return altSplicJunction;
    }

    private static final int MAX_HOMOLOGY_LENGTH = 5;

    private void checkJunctionHomology(final int[] spliceJunction, final Map<RegionReadData,RegionMatchType> readRegions)
    {
        if(readRegions.values().stream().anyMatch(x -> x == EXON_BOUNDARY)) // only consider if neither side matches a known junction
            return;

        // check for a junction position (or both) which can be moved by homology to match a known splice site
        final String[] junctionBases = new String[] {"", ""}; // 5 bases either side of the junction
        int[] nearestExonBoundary = {0, 0}; // +ve if SJ is beyond the boundary

        for(int se = SE_START; se <= SE_END; ++se)
        {
            junctionBases[se] = mConfig.RefGenome.getBaseString(
                    mGenes.chromosome(), spliceJunction[se] - MAX_HOMOLOGY_LENGTH, spliceJunction[se] + MAX_HOMOLOGY_LENGTH);

            int nearestBoundary = MAX_HOMOLOGY_LENGTH + 1;

            for(RegionReadData region : readRegions.keySet())
            {
                int distance = abs(spliceJunction[se] - region.start()) < abs(spliceJunction[se] - region.end()) ?
                        spliceJunction[se] - region.start() : spliceJunction[se] - region.end();

                if(abs(distance) < abs(nearestBoundary))
                    nearestBoundary = distance;
            }

            nearestExonBoundary[se] = nearestBoundary;
        }

        if(abs(nearestExonBoundary[SE_START]) > MAX_HOMOLOGY_LENGTH && abs(nearestExonBoundary[SE_END]) > MAX_HOMOLOGY_LENGTH)
            return;

        // check the closer of the 2 ends
        int minDistance = min(nearestExonBoundary[SE_START], nearestExonBoundary[SE_END]);
        int juncBaseIndex = MAX_HOMOLOGY_LENGTH;

        if(minDistance > 0)
        {
            // SJ is past the nearest exon boundary, so compare the preceding bases
            final String startBases = junctionBases[SE_START].substring(juncBaseIndex - abs(minDistance), juncBaseIndex);
            final String endBases = junctionBases[SE_START].substring(juncBaseIndex - abs(minDistance), juncBaseIndex);

            if(!startBases.equals(endBases))
                return;
        }
        else
        {
            // try to shift the junction forwards - requiring the bases to the right at the start to match at each end

            // example of the SJ being 2 bases back from the exon boundary
            // 01234567
            //      S E

            final String startBases = junctionBases[SE_START].substring(juncBaseIndex, juncBaseIndex + abs(minDistance));
            final String endBases = junctionBases[SE_START].substring(juncBaseIndex, juncBaseIndex + abs(minDistance));

            if(!startBases.equals(endBases))
                return;
        }

        // shift the bases
        spliceJunction[SE_START] -= minDistance;
        spliceJunction[SE_END] -= minDistance;
    }

    public AltSpliceJunction createFromReads(final ReadRecord read1, final ReadRecord read2, final List<Integer> relatedTransIds)
    {
        // related transcripts will any of those where either read covers one or more of its exons
        int[] spliceJunction = new int[SE_PAIR];

        List<RegionReadData> sjStartRegions = Lists.newArrayList(); // transcript regions with an exon matching the start of the alt SJ
        List<RegionReadData> sjEndRegions = Lists.newArrayList();
        List<RegionReadData> candidateRegions = Lists.newArrayList();

        for(int i = 0; i <= 1; ++i)
        {
            final ReadRecord read = (i == 0) ? read1 : read2;

            // take the longer of the 2 soft-clippings
            ClippedSide clippedSide = ClippedSide.from(read.Cigar, read.isSoftClipped(SE_START), read.isSoftClipped(SE_END));

            boolean useLeft = clippedSide != null && clippedSide.isLeft();

            if(useLeft)
                spliceJunction[SE_START] = read.getCoordsBoundary(SE_START);
            else
                spliceJunction[SE_END] = read.getCoordsBoundary(SE_END);

            for(RegionReadData region : read.getMappedRegions().keySet())
            {
                if(!candidateRegions.contains(region))
                    candidateRegions.add(region);
            }
        }

        if(spliceJunction[SE_START] == 0 || spliceJunction[SE_END] == 0)
            return null;

        AltSpliceJunction altSplicJunction;

        if(spliceJunction[SE_START] < spliceJunction[SE_END])
        {
            // flip the SJ around since a DUP/circular has opposite orientations to a standard alt-SJ DL
            int[] flippedSpliceJunction = { spliceJunction[SE_END], spliceJunction[SE_START] };
            final AltSpliceJunctionContext[] flippedRegionContexts = { AltSpliceJunctionContext.UNKNOWN, AltSpliceJunctionContext.UNKNOWN };

            classifyRegions(read1, flippedSpliceJunction, sjEndRegions, sjStartRegions, flippedRegionContexts);
            classifyRegions(read2, flippedSpliceJunction, sjEndRegions, sjStartRegions, flippedRegionContexts);

            final AltSpliceJunctionContext[] regionContexts = { flippedRegionContexts[SE_END], flippedRegionContexts[SE_START] };

            altSplicJunction = new AltSpliceJunction(
                    mGenes.chromosome(), spliceJunction, CIRCULAR, read1.Id, regionContexts, sjStartRegions, sjEndRegions);
        }
        else
        {
            // pair of reads with a DEL-type split but using supplementary data - treat normally
            int pos = spliceJunction[SE_END];
            spliceJunction[SE_END] = spliceJunction[SE_START];
            spliceJunction[SE_START] = pos;

            final AltSpliceJunctionContext[] regionContexts = { AltSpliceJunctionContext.UNKNOWN, AltSpliceJunctionContext.UNKNOWN};

            classifyRegions(read1, spliceJunction, sjEndRegions, sjStartRegions, regionContexts);
            classifyRegions(read2, spliceJunction, sjEndRegions, sjStartRegions, regionContexts);

            AltSpliceJunctionType sjType = classifySpliceJunction(relatedTransIds, sjStartRegions, sjEndRegions, regionContexts);

            altSplicJunction = new AltSpliceJunction(
                    mGenes.chromosome(), spliceJunction, sjType, read1.Id, regionContexts, sjStartRegions, sjEndRegions);
        }

        altSplicJunction.setCandidateTranscripts(candidateRegions);

        return altSplicJunction;
    }

    private void classifyRegions(
            final ReadRecord read, final int[] spliceJunction,
            final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions, final AltSpliceJunctionContext[] regionContexts)
    {
        // collect up all exon regions matching the observed novel splice junction
        final List<Integer> sjMatchedTransIds = Lists.newArrayList();

        for (Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
        {
            final RegionReadData region = entry.getKey();
            RegionMatchType matchType = entry.getValue();

            if(matchType == RegionMatchType.NONE)
                continue;

            if (region.end() == spliceJunction[SE_START])
            {
                // the region cannot only be from a transcript ending here
                if(region.getPostRegions().isEmpty())
                    continue;

                regionContexts[SE_START] = SPLICE_JUNC;
                sjStartRegions.add(region);
                region.getTransExonRefs().forEach(x -> sjMatchedTransIds.add(x.TransId));
            }

            if (region.start() == spliceJunction[SE_END])
            {
                if(region.getPreRegions().isEmpty())
                    continue;

                regionContexts[SE_END] = SPLICE_JUNC;
                sjEndRegions.add(region);
                region.getTransExonRefs().forEach(x -> sjMatchedTransIds.add(x.TransId));
            }
        }

        // now check for exonic
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(regionContexts[se] == SPLICE_JUNC)
                continue;

            for (Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
            {
                final RegionReadData region = entry.getKey();
                RegionMatchType matchType = entry.getValue();

                if (matchType == RegionMatchType.NONE)
                    continue;

                if(!sjMatchedTransIds.isEmpty() && !region.getTransExonRefs().stream().anyMatch(x -> sjMatchedTransIds.contains(x.TransId)))
                    continue;

                if(positionWithin(spliceJunction[se], region.start(), region.end()))
                {
                    regionContexts[se] = EXONIC;
                }
            }
        }

        // and any remaining to assign as intronic
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(regionContexts[se] == AltSpliceJunctionContext.UNKNOWN)
                regionContexts[se] = AltSpliceJunctionContext.INTRONIC;
        }
    }

    private AltSpliceJunctionType classifySpliceJunction(
            final List<Integer> transIds, final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions,
            final AltSpliceJunctionContext[] regionContexts)
    {
        /* classify the overall type of novel splice junction
            - Exon_skipping - both sides match known SJ
            - Novel 5’ or 3' splice site (one side matches known SJ
            - Novel exon - phased Novel 5’ and 3’ splice site from 2 reads in the same fragment
            - Novel intron - both sides exonic - may also be a transcribed somatic DEL
            - Intronic - both sides Intronic  - may also be a transcribed somatic DEL
        */

        if(!sjStartRegions.isEmpty() || !sjEndRegions.isEmpty())
        {
            boolean hasStartMatch = false;
            boolean hasEndMatch = false;
            int startTransId = -1;
            int endTransId = -1;

            for(Integer transId : transIds)
            {
                // check for skipped exons indicated by the same transcript matching the start and end of the alt SJ, but skipping an exon
                boolean matchesStart = sjStartRegions.stream().anyMatch(x -> x.hasTransId(transId));
                boolean matchesEnd = sjEndRegions.stream().anyMatch(x -> x.hasTransId(transId));

                if(matchesStart && matchesEnd)
                    return SKIPPED_EXONS;

                if(matchesStart)
                {
                    hasStartMatch = true;
                    startTransId = transId;
                }
                else if(matchesEnd)
                {
                    hasEndMatch = true;
                    endTransId = transId;
                }
            }

            if(hasStartMatch && hasEndMatch) // 2 different transcripts match the SJs
                return MIXED_TRANS;

            if(hasStartMatch)
            {
                boolean forwardStrand = mGenes.getStrand(startTransId) == 1;
                return forwardStrand ? NOVEL_3_PRIME : NOVEL_5_PRIME;
            }
            else if(hasEndMatch)
            {
                boolean forwardStrand = mGenes.getStrand(endTransId) == 1;
                return forwardStrand ? NOVEL_5_PRIME : NOVEL_3_PRIME;
            }
        }

        if(regionContexts[SE_START] == AltSpliceJunctionContext.INTRONIC && regionContexts[SE_END] == AltSpliceJunctionContext.INTRONIC)
            return INTRONIC;

        if(regionContexts[SE_START] == EXONIC && regionContexts[SE_END] == EXONIC)
            return NOVEL_INTRON;

        return EXON_INTRON;
    }

    public void checkNovelExon(AltSpliceJunction firstAltSJ, AltSpliceJunction secondAltSJ)
    {
        if(!(firstAltSJ.type() == NOVEL_3_PRIME && secondAltSJ.type() == NOVEL_5_PRIME)
        && !(firstAltSJ.type() == NOVEL_5_PRIME && secondAltSJ.type() == NOVEL_3_PRIME))
        {
            return;
        }

        final List<RegionReadData> regions1 = Lists.newArrayList(firstAltSJ.getSjStartRegions());
        regions1.addAll(firstAltSJ.getSjEndRegions());

        final List<RegionReadData> regions2 = Lists.newArrayList(secondAltSJ.getSjStartRegions());
        regions2.addAll(secondAltSJ.getSjEndRegions());

        List<Integer> commonTranscripts = Lists.newArrayList();

        for(final RegionReadData region1 : regions1)
        {
            List<Integer> transIds1 = region1.getTransExonRefs().stream().map(x -> x.TransId).collect(Collectors.toList());

            for(Integer transId1 : transIds1)
            {
                for(final RegionReadData region2 : regions2)
                {
                    List<Integer> transIds2 = region2.getTransExonRefs().stream().map(x -> x.TransId).collect(Collectors.toList());

                    if (transIds2.contains(transId1))
                    {
                        if(!commonTranscripts.contains(transId1))
                            commonTranscripts.add(transId1);

                        break;
                    }
                }
            }
        }

        if(commonTranscripts.isEmpty())
            return;

        firstAltSJ.overrideType(NOVEL_EXON);
        secondAltSJ.overrideType(NOVEL_EXON);

        // remove any transcripts not supported by both reads
        firstAltSJ.cullNonMatchedTranscripts(commonTranscripts);
        secondAltSJ.cullNonMatchedTranscripts(commonTranscripts);
    }

    private AltSpliceJunction registerAltSpliceJunction(
            final List<GeneReadData> candidateGenes, final ReadRecord read, final List<Integer> regionTranscripts)
    {
        AltSpliceJunction altSpliceJunc = createFromRead(read, regionTranscripts);

        AltSpliceJunction existingSpliceJunc = mAltSpliceJunctions.stream()
                .filter(x -> x.matches(altSpliceJunc)).findFirst().orElse(null);

        if(existingSpliceJunc != null)
        {
            existingSpliceJunc.addFragmentCount();
            return existingSpliceJunc;
        }

        // extra check that the SJ doesn't match any transcript
        if(junctionMatchesKnownSpliceJunction(candidateGenes, altSpliceJunc.SpliceJunction, regionTranscripts))
            return null;

        altSpliceJunc.addFragmentCount();

        mAltSpliceJunctions.add(altSpliceJunc);
        return altSpliceJunc;
    }

    private AltSpliceJunction registerAltSpliceJunction(
            final List<GeneReadData> candidateGenes, final ReadRecord read1, final ReadRecord read2, final List<Integer> regionTranscripts)
    {
        AltSpliceJunction altSpliceJunc = createFromReads(read1, read2, regionTranscripts);

        if(altSpliceJunc == null)
            return null;

        AltSpliceJunction existingSpliceJunc = mAltSpliceJunctions.stream()
                .filter(x -> x.matches(altSpliceJunc)).findFirst().orElse(null);

        if(existingSpliceJunc != null)
        {
            existingSpliceJunc.addFragmentCount();
            return existingSpliceJunc;
        }

        if(junctionMatchesKnownSpliceJunction(candidateGenes, altSpliceJunc.SpliceJunction, regionTranscripts))
            return null;

        altSpliceJunc.addFragmentCount();

        mAltSpliceJunctions.add(altSpliceJunc);
        return altSpliceJunc;
    }

    public void setPositionDepth(final BaseDepth baseDepth)
    {
        for(AltSpliceJunction altSJ : mAltSpliceJunctions)
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                int depth = baseDepth.depthAtBase(altSJ.SpliceJunction[se]);
                altSJ.addPositionCount(se, depth);
            }
        }
    }

    public void prioritiseGenes()
    {
        mAltSpliceJunctions.forEach(x -> x.setBaseContext(mConfig.RefGenome, mGenes.chromosome()));

        for(AltSpliceJunction altSJ : mAltSpliceJunctions)
        {
            final List<Integer> transIds = altSJ.candidateTransIds();

            GeneReadData topGene = null;
            int topMatch = 0;

            int spliceStrand = altSJ.getKnownSpliceBaseStrand();

            List<GeneReadData> candidateGenes = Lists.newArrayList();

            for(final GeneReadData gene : mGenes.genes())
            {
                if(gene.getExonRegions().stream().anyMatch(x -> altSJ.getSjStartRegions().contains(x) || altSJ.getSjEndRegions().contains(x)))
                {
                    candidateGenes.add(gene);
                }
            }

            if(candidateGenes.isEmpty())
                candidateGenes = mGenes.genes();

            for(final GeneReadData gene : candidateGenes)
            {
                if(spliceStrand != 0 && gene.GeneData.Strand != spliceStrand)
                    continue;

                int transMatched = (int)gene.getTranscripts().stream().filter(x -> transIds.contains(x.TransId)).count();

                if(transMatched > topMatch)
                {
                    topMatch = transMatched;
                    topGene = gene;
                }
            }

            if(topGene != null)
            {
                altSJ.setGeneId(topGene.GeneData.GeneId);
            }
            else
            {
                altSJ.setGeneId(mGenes.genes().get(0).GeneData.GeneId);
            }
        }
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile(ALT_SJ_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(AltSpliceJunction.csvHeader());
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create alt splice junction writer: {}", e.toString());
            return null;
        }
    }

    public void writeAltSpliceJunctions()
    {
        if(mWriter != null)
        {
            for(final AltSpliceJunction altSJ : mAltSpliceJunctions)
            {
                final GeneReadData gene = mGenes.genes().stream().filter(x -> x.GeneData.GeneId.equals(altSJ.getGeneId())).findFirst().orElse(null);
                altSJ.calcSummaryData(gene);
            }

            writeAltSpliceJunctions(mWriter, mAltSpliceJunctions, mGenes);
        }
    }

    private synchronized static void writeAltSpliceJunctions(
            final BufferedWriter writer, final List<AltSpliceJunction> altSpliceJunctions, final GeneCollection geneCollection)
    {
        try
        {
            for(final AltSpliceJunction altSJ : altSpliceJunctions)
            {
                final GeneReadData gene = geneCollection.genes().stream()
                        .filter(x -> x.GeneData.GeneId.equals(altSJ.getGeneId())).findFirst().orElse(null);

                writer.write(altSJ.toCsv(gene.GeneData));
                writer.newLine();
            }

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt splice junction file: {}", e.toString());
        }
    }

}
