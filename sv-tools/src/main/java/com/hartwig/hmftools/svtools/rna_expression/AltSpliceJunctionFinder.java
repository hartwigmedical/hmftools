package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_EXONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_SJ;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.MIXED_TRANS;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.generateMappedCoords;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransId;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.deriveCommonRegions;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionWithin;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class AltSpliceJunctionFinder
{
    private final RnaExpConfig mConfig;
    private final SamReader mSamReader;

    private final List<AltSpliceJunction> mAltSpliceJunctions;
    private final BufferedWriter mWriter;

    private GeneReadData mGene;
    private long mCurrentSjPosition;

    public AltSpliceJunctionFinder(final RnaExpConfig config, final SamReader samReader, final BufferedWriter writer)
    {
        mConfig = config;
        mSamReader = samReader;
        mAltSpliceJunctions = Lists.newArrayList();
        mWriter = writer;
        mCurrentSjPosition = 0;
        mGene = null;
    }

    public List<AltSpliceJunction> getAltSpliceJunctions() { return mAltSpliceJunctions; }

    public void setGeneData(final GeneReadData gene)
    {
        mGene = gene;
        mAltSpliceJunctions.clear();
        mFragmentTracker.clear();
    }

    public void evaluateFragmentReads(final ReadRecord read1, final ReadRecord read2, final List<Integer> relatedTransIds)
    {
        if(mGene == null)
            return;

        // for now exclude SJs outside known transcripts
        if(read1.PosStart < mGene.GeneData.GeneStart || read2.PosStart < mGene.GeneData.GeneStart
        || read1.PosEnd > mGene.GeneData.GeneEnd || read2.PosEnd > mGene.GeneData.GeneEnd)
        {
            return;
        }

        AltSpliceJunction firstAltSJ = null;

        if(AltSpliceJunctionFinder.isCandidate(read1))
        {
            firstAltSJ = registerAltSpliceJunction(read1, relatedTransIds);
        }

        AltSpliceJunction secondAltSJ = null;

        if(AltSpliceJunctionFinder.isCandidate(read2))
        {
            // avoid double-counting overlapping reads
            if(firstAltSJ != null && positionsOverlap(read1.PosStart, read1.PosEnd, read2.PosStart, read2.PosEnd))
                return;

            secondAltSJ = registerAltSpliceJunction(read2, relatedTransIds);
        }

        if(firstAltSJ != null && secondAltSJ != null)
        {
            checkNovelExon(firstAltSJ, secondAltSJ);
        }
    }

    private static boolean isCandidate(final ReadRecord read)
    {
        if(!read.Cigar.containsOperator(CigarOperator.N))
            return false;

        if(read.getTranscriptClassifications().values().contains(SPLICE_JUNCTION))
            return false;

        if(read.getMappedRegionCoords().size() == 1)
            return false;

        return true;
    }

    public boolean junctionMatchesGene(final long[] spliceJunction, final List<Integer> transIds)
    {
        for(TranscriptData transData : mGene.getTranscripts())
        {
            if(!transIds.contains(transData.TransId))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                if(transData.exons().get(i).ExonEnd == spliceJunction[SE_START] && transData.exons().get(i+1).ExonStart == spliceJunction[SE_END])
                    return true;
            }
        }

        // also check any overlapping gene's exonic regions
        if(mGene.getOtherGeneExonicRegions().stream()
                .anyMatch(x -> positionsOverlap(spliceJunction[SE_START], spliceJunction[SE_END], x[SE_START], x[SE_END])))
        {
            return true;
        }

        return false;
    }

    public AltSpliceJunction createFromRead(final ReadRecord read, final List<Integer> relatedTransIds)
    {
        // related transcripts will any of those where either read covers one or more of its exons
        long[] spliceJunction = new long[SE_PAIR];

        // find the novel splice junction, and all associated transcripts
        final List<long[]> mappedCoords = read.getMappedRegionCoords();

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

        List<RegionReadData> sjStartRegions = Lists.newArrayList(); // transcript regions with an exon matching the start of the alt SJ
        List<RegionReadData> sjEndRegions = Lists.newArrayList();
        String[] regionContexts = {"", ""};

        classifyRegions(read, spliceJunction, sjStartRegions, sjEndRegions, regionContexts);

        AltSpliceJunctionType sjType = classifySpliceJunction(relatedTransIds, sjStartRegions, sjEndRegions, regionContexts);

        AltSpliceJunction altSplicJunction = new AltSpliceJunction(mGene, spliceJunction, sjType, regionContexts);

        altSplicJunction.StartRegions.addAll(sjStartRegions);
        altSplicJunction.EndRegions.addAll(sjEndRegions);

        altSplicJunction.setCandidateTranscripts(read.getMappedRegions().keySet().stream().collect(Collectors.toList()));

        return altSplicJunction;
    }

    private void classifyRegions(
            final ReadRecord read, final long[] spliceJunction,
            final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions, String[] regionContexts)
    {
        // collect up all exon regions matching the observed novel splice junction
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

                regionContexts[SE_START] = CONTEXT_SJ;
                sjStartRegions.add(region);
            }
            else if(positionWithin(spliceJunction[SE_START], region.start(), region.end()))
            {
                if(regionContexts[SE_START] != CONTEXT_SJ)
                    regionContexts[SE_START] = CONTEXT_EXONIC;
            }

            if (region.start() == spliceJunction[SE_END])
            {
                if(region.getPreRegions().isEmpty())
                    continue;

                regionContexts[SE_END] = CONTEXT_SJ;
                sjEndRegions.add(region);
            }
            else if(positionWithin(spliceJunction[SE_END], region.start(), region.end()))
            {
                if(regionContexts[SE_END] != CONTEXT_SJ)
                    regionContexts[SE_END] = CONTEXT_EXONIC;
            }
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(regionContexts[se] == "")
                regionContexts[se] = CONTEXT_INTRONIC;
        }
    }

    private AltSpliceJunctionType classifySpliceJunction(
            final List<Integer> transIds, final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions,
            final String[] regionContexts)
    {
        /* classify the overall type of novel splice junction
            - Exon_skipping - both sides match known SJ
            - Novel 5’ or 3' splice site (one side matches known SJ
            - Novel exon - phased Novel 5’ and 3’ splice site from 2 reads in the same fragment
            - Novel intron - both sides exonic - may also be a transcribed somatic DEL
            - Intronic - both sides Intronic  - may also be a transcribed somatic DEL
        */

        if (!sjStartRegions.isEmpty() || !sjEndRegions.isEmpty())
        {
            boolean hasStartMatch = false;
            boolean hasEndMatch = false;

            for (Integer transId : transIds)
            {
                // check for skipped exons indicated by the same transcript matching the start and end of the alt SJ, but skipping an exon
                boolean matchesStart = sjStartRegions.stream().anyMatch(x -> x.hasTransId(transId));
                boolean matchesEnd = sjEndRegions.stream().anyMatch(x -> x.hasTransId(transId));

                if (matchesStart && matchesEnd)
                    return SKIPPED_EXONS;

                hasStartMatch |= matchesStart;
                hasEndMatch |= matchesEnd;
            }

            if(hasStartMatch && hasEndMatch) // 2 different transcripts match the SJs
                return MIXED_TRANS;
            else if(hasStartMatch)
                return mGene.GeneData.forwardStrand() ? NOVEL_3_PRIME : NOVEL_5_PRIME;
            else if(hasEndMatch)
                return mGene.GeneData.forwardStrand() ? NOVEL_5_PRIME : NOVEL_3_PRIME;
        }

        if(regionContexts[SE_START] == CONTEXT_INTRONIC && regionContexts[SE_END] == CONTEXT_INTRONIC)
            return INTRONIC;

        if(regionContexts[SE_START] == CONTEXT_EXONIC && regionContexts[SE_END] == CONTEXT_EXONIC)
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

        final List<RegionReadData> regions1 = Lists.newArrayList(firstAltSJ.StartRegions);
        regions1.addAll(firstAltSJ.EndRegions);

        final List<RegionReadData> regions2 = Lists.newArrayList(secondAltSJ.StartRegions);
        regions2.addAll(secondAltSJ.EndRegions);

        List<Integer> commonTranscripts = Lists.newArrayList();

        for(final RegionReadData region1 : regions1)
        {
            List<Integer> transIds1 = region1.getRefRegions().stream().map(x -> extractTransId(x)).collect(Collectors.toList());

            for(Integer transId1 : transIds1)
            {
                for(final RegionReadData region2 : regions2)
                {
                    List<Integer> transIds2 = region2.getRefRegions().stream().map(x -> extractTransId(x)).collect(Collectors.toList());

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

    private AltSpliceJunction registerAltSpliceJunction(final ReadRecord read, final List<Integer> regionTranscripts)
    {
        if(mGene == null)
            return null;

        AltSpliceJunction altSpliceJunc = createFromRead(read, regionTranscripts);

        AltSpliceJunction existingSpliceJunc = mAltSpliceJunctions.stream()
                .filter(x -> x.matches(altSpliceJunc)).findFirst().orElse(null);

        if(existingSpliceJunc != null)
        {
            existingSpliceJunc.addFragmentCount();
            return existingSpliceJunc;
        }

        // extra check that the SJ doesn't match any transcript
        if(junctionMatchesGene(altSpliceJunc.SpliceJunction, regionTranscripts))
            return null;

        altSpliceJunc.addFragmentCount();

        mAltSpliceJunctions.add(altSpliceJunc);
        return altSpliceJunc;
    }

    // read depth count state
    private PerformanceCounter mReadDepthPerf = new PerformanceCounter("AltSJ ReadDepth");
    private FragmentTracker mFragmentTracker = new FragmentTracker();
    private int mTotalReadsProcessed = 0;

    public void recordDepthCounts(int totalReadCount)
    {
        if(mAltSpliceJunctions.isEmpty())
            return;

        BamSlicer slicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        QueryInterval[] queryInterval = new QueryInterval[1];
        int chrSeqIndex = mSamReader.getFileHeader().getSequenceIndex(mGene.GeneData.Chromosome);

        mFragmentTracker.clear();
        mTotalReadsProcessed = 0;

        if(RE_LOGGER.isDebugEnabled() && (mAltSpliceJunctions.size() >= 50 || totalReadCount > 100000))
        {
            long totalSJRange =
                    mAltSpliceJunctions.stream().mapToLong(x -> x.SpliceJunction[SE_END] - x.SpliceJunction[SE_START]).sum();
            double avgSJLength = totalSJRange / (double) mAltSpliceJunctions.size();
            long geneLength = mGene.GeneData.GeneEnd - mGene.GeneData.GeneStart;
            double expReadsPerSJ = avgSJLength / (double) geneLength * totalReadCount;

            RE_LOGGER.debug(String.format("gene(%s) length(%d) totalReads(%d) altSJs(count=%d totalLen=%d avgLen=%.0f) expReads(%.0f)",
                    mGene.name(), geneLength, totalReadCount, mAltSpliceJunctions.size(), totalSJRange, avgSJLength, expReadsPerSJ));
        }

        // String geneDetails = String.format("%s,%d,%d", mGene.name(), totalReadCount, mAltSpliceJunctions.size());
        // mReadDepthPerf.start(geneDetails);
        mReadDepthPerf.start();

        // for genes of super high depth, take all reads for the gene together and then test each against each SJ
        long altSJMinPos = mAltSpliceJunctions.stream().mapToLong(x -> x.SpliceJunction[SE_START]).min().orElse(mGene.GeneData.GeneStart);
        long altSJMaxPos = mAltSpliceJunctions.stream().mapToLong(x -> x.SpliceJunction[SE_END]).max().orElse(mGene.GeneData.GeneEnd);

        queryInterval[0] = new QueryInterval(chrSeqIndex, (int)altSJMinPos, (int)altSJMaxPos);

        slicer.slice(mSamReader, queryInterval, this::setPositionDepthFromRead);

        mReadDepthPerf.stop();

        // RE_LOGGER.debug("gene({}) read depth counts altSJs({}) readsProcessed({}) unmatchedFragReads({})",
        //        mGene.name(), mAltSpliceJunctions.size(), mTotalReadsProcessed, mFragmentTracker.readsCount());
    }

    public void logStats()
    {
        mReadDepthPerf.logStats();
    }

    private void setPositionDepthFromRead(@NotNull final SAMRecord record)
    {
        final List<long[]> readCoords = generateMappedCoords(record.getCigar(), record.getStart());
        final List<long[]> otherReadCoords = (List<long[]>)mFragmentTracker.checkRead(record.getReadName(), readCoords);

        if(otherReadCoords == null)
            return;

        ++mTotalReadsProcessed;

        final List<long[]> commonMappings = deriveCommonRegions(readCoords, otherReadCoords);
        long readMinPos = commonMappings.get(0)[SE_START];
        long readMaxPos = commonMappings.get(commonMappings.size() - 1)[SE_END];

        for(AltSpliceJunction altSJ : mAltSpliceJunctions)
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                int position = (int) altSJ.SpliceJunction[se];

                if(!positionWithin(position, readMinPos, readMaxPos))
                    continue;

                if (commonMappings.stream().anyMatch(x -> positionWithin(position, x[SE_START], x[SE_END])))
                {
                    altSJ.addPositionCount(se);
                }
            }
        }
    }

    public static BufferedWriter createAltSpliceJunctionWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("alt_splice_junc.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,Chromosome,Strand,SjStart,SjEnd,FragCount,StartDepth,EndDepth");
            writer.write(",Type,StartContext,EndContext,NearestStartExon,NearestEndExon");
            writer.write(",StartBases,EndBases,StartTrans,EndTrans");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            RE_LOGGER.error("failed to create at splice junction writer: {}", e.toString());
            return null;
        }
    }

    public void writeAltSpliceJunctions()
    {
        if(mWriter != null)
            writeAltSpliceJunctions(mWriter, mAltSpliceJunctions, mConfig);
    }

    private synchronized static void writeAltSpliceJunctions(
            final BufferedWriter writer, final List<AltSpliceJunction> altSpliceJunctions, final RnaExpConfig config)
    {
        try
        {
            for(final AltSpliceJunction altSJ : altSpliceJunctions)
            {
                writer.write(String.format("%s,%s,%s,%d",
                        altSJ.Gene.GeneData.GeneId, altSJ.Gene.GeneData.GeneName,
                        altSJ.Gene.GeneData.Chromosome, altSJ.Gene.GeneData.Strand));

                writer.write(String.format(",%d,%d,%d,%d,%d",
                        altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END], altSJ.getFragmentCount(),
                        altSJ.getPositionCount(SE_START), altSJ.getPositionCount(SE_END)));

                final String startTrans = altSJ.StartRegions.isEmpty() ? "NONE" : altSJ.startTranscriptNames();
                final String endTrans = altSJ.EndRegions.isEmpty() ? "NONE" : altSJ.endTranscriptNames();

                writer.write(String.format(",%s,%s,%s,%d,%d,%s,%s,%s,%s",
                        altSJ.type(), altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                        altSJ.calcNearestExonBoundary(SE_START), altSJ.calcNearestExonBoundary(SE_END),
                        altSJ.getBaseContext(config.RefFastaSeqFile, SE_START), altSJ.getBaseContext(config.RefFastaSeqFile, SE_END),
                        startTrans, endTrans));

                writer.newLine();
            }

        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write alt splice junction file: {}", e.toString());
        }
    }



}
