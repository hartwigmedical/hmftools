package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.biomodalcollapse.AlignmentStats.getAlignmentStats;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapse.writeStatLine;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.QualCappingOption.CAP_BY_FIRST;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.collapseAlignedSeqExact;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.collapseAlignedSeqModC;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.consensusReadForStatOutput;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.fastqToSeq;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.getCigar;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.getCutPoint;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.getExactScore;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.getLeftElements;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.getRightElements;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.qualityTrim;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.reverseComplementSeq;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.sanatizeQualString;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.seqToFastq;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.zip;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.EXTEND_GAP_PENALTY;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.FORWARD_HAIRPIN;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.INS_BASE;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.MIN_RESOLVED_READ_LENGTH;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.MISSING_BASE;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.MODC_BASE;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.OPEN_GAP_PENALTY;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.REVERSE_HAIRPIN;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.STAT_DELIMITER;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.STAT_HEADERS;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.HairpinInfo.findHairpin;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.ReverseComplementMatchInfo.findBestReverseComplementMatch;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.BufferedWriter;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqRecord;

public class BiomodalCollapseWorker extends Thread
{
    private static final BaseQualPair INS_BASE_QUAL_QUAL = new BaseQualPair(INS_BASE, 0);

    public final int ThreadId;

    private final SynchronizedPairedFastqReader mFastqPairReader;
    private final BufferedWriter mResolvedFastqWriter;
    private final BufferedWriter mDebugStatsWriter;
    private final Map<String, FastqRecord> mRefResolvesFastqMap;
    private final BiomodalCollapseStats mStats;
    private final SynchronizedPairedFastqWriter mBadFastqPairWriter;

    private final NeedlemanWunschAligner<BaseQualPair> mAligner;

    public BiomodalCollapseWorker(int threadId, final SynchronizedPairedFastqReader fastqPairReader,
            final BufferedWriter resolvedFastqWriter, @Nullable final BufferedWriter debugStatsWriter,
            final Map<String, FastqRecord> refResolvesFastqMap, final BiomodalCollapseStats stats,
            @Nullable final SynchronizedPairedFastqWriter badFastqPairWriter)
    {
        ThreadId = threadId;
        mFastqPairReader = fastqPairReader;
        mResolvedFastqWriter = resolvedFastqWriter;
        mDebugStatsWriter = debugStatsWriter;
        mStats = stats;
        mRefResolvesFastqMap = refResolvesFastqMap;
        mBadFastqPairWriter = badFastqPairWriter;

        mAligner = new NeedlemanWunschAligner<>();
    }

    @Override
    public void run()
    {
        Pair<FastqRecord, FastqRecord> fastqPair;
        while((fastqPair = mFastqPairReader.getNext()) != null)
        {
            FastqRecord resolvedFastq = processFastqPair(fastqPair.getLeft(), fastqPair.getRight());

            int processedFastqPairs = mStats.ProcessedFastqPairCount.incrementAndGet();
            if(processedFastqPairs % 10_000 == 0)
            {
                BT_LOGGER.info(format("%d fastq pairs have been processed", processedFastqPairs));
            }

            if(resolvedFastq != null)
            {
                mStats.WrittenConsensusFastqCount.getAndIncrement();
                BiomodalCollapse.writeResolvedFastqRecord(mResolvedFastqWriter, resolvedFastq);
            }
        }
    }

    @Nullable
    private FastqRecord processFastqPair(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        int trimmedLength = min(fastq1.getReadLength(), fastq2.getReadLength());

        final List<BaseQualPair> seq1 = fastqToSeq(fastq1, 0, trimmedLength - 1);
        final List<BaseQualPair> seq2 = fastqToSeq(fastq2, 0, trimmedLength - 1);

        final List<BaseQualPair> seq1RC = reverseComplementSeq(seq1);
        final List<BaseQualPair> seq2RC = reverseComplementSeq(seq2);

        ReverseComplementMatchInfo rcMatch = findBestReverseComplementMatch(seq1, seq2RC);

        HairpinInfo hairpin1 = findHairpin(seq1, FORWARD_HAIRPIN);
        HairpinInfo hairpin2 = findHairpin(seq2, REVERSE_HAIRPIN);

        int cutPoint = getCutPoint(trimmedLength, hairpin1, hairpin2, rcMatch);
        if(cutPoint <= 0)
        {
            if(mDebugStatsWriter != null)
            {
                writeStats(fastq1, fastq2, hairpin1, hairpin2, rcMatch, seq1, seq2, seq1RC, seq2RC, null, null,
                        null, null, 0, null, null, 0,
                        null, 0, 0);
            }

            mStats.BadWrittenFastqPairCount.getAndIncrement();
            if(mBadFastqPairWriter != null)
            {
                mBadFastqPairWriter.write(fastq1, fastq2);
            }

            return null;
        }

        // forward aligned consensus
        List<BaseQualPair> read1 = seq1.subList(0, cutPoint);
        List<BaseQualPair> read2 = seq2.subList(0, cutPoint);
        List<Pair<BaseQualPair, BaseQualPair>> forwardAlignment =
                mAligner.approxAlign(read1, read2, BiomodalCollapseUtil::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, true, true);
        List<BaseQualPair> forwardConsensus = collapseAlignedSeqModC(forwardAlignment);

        // reverse aligned consensus
        int rcEndIndex = -1;
        List<Pair<BaseQualPair, BaseQualPair>> reverseAlignment = null;
        List<BaseQualPair> reverseConsensus = null;
        if(rcMatch != null)
        {
            int rcStartIndex = -rcMatch.Read1Shift;
            int rcLength = min(trimmedLength, cutPoint - rcStartIndex);
            rcEndIndex = rcStartIndex + rcLength - 1;
            if(rcLength > 0)
            {
                List<BaseQualPair> read1RC = seq1RC.subList(0, rcLength);
                List<BaseQualPair> read2RC = seq2RC.subList(0, rcLength);
                reverseAlignment =
                        mAligner.approxAlign(read2RC, read1RC, BiomodalCollapseUtil::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, true, false, false);
                reverseConsensus = collapseAlignedSeqModC(reverseAlignment);
            }
        }

        // final consensus
        int reverseConsensusTrimCount = 0;
        List<Pair<BaseQualPair, BaseQualPair>> finalAlignment = null;
        if(reverseAlignment == null)
        {
            finalAlignment = zip(forwardConsensus, Collections.emptyList());
        }
        else
        {
            // this is expensive to fully align, and indels are rare, so we naively align
            Collections.reverse(forwardConsensus);
            Collections.reverse(reverseConsensus);
            finalAlignment = zip(forwardConsensus, reverseConsensus);
            Collections.reverse(forwardConsensus);
            Collections.reverse(reverseConsensus);
            Collections.reverse(finalAlignment);

            for(int i = 0; i < finalAlignment.size(); i++)
            {
                if(finalAlignment.get(i).getLeft() != null)
                {
                    break;
                }

                reverseConsensusTrimCount++;
            }

            finalAlignment = finalAlignment.subList(reverseConsensusTrimCount, finalAlignment.size());
        }

        List<BaseQualPair> finalConsensus = collapseAlignedSeqExact(finalAlignment, CAP_BY_FIRST);

        // clipping final consensus
        Pair<Integer, Integer> clippingBounds = qualityTrim(finalConsensus);
        if(clippingBounds != null && clippingBounds.getRight() - clippingBounds.getLeft() + 1 < MIN_RESOLVED_READ_LENGTH)
        {
            clippingBounds = null;
        }

        int prefixTrimCount = clippingBounds == null ? finalConsensus.size() : clippingBounds.getLeft();
        int suffixTrimCount = clippingBounds == null ? 0 : finalConsensus.size() - 1 - clippingBounds.getRight();

        FastqRecord finalFastq = null;
        List<BaseQualPair> clippedConsensus = null;
        if(clippingBounds != null)
        {
            clippedConsensus = finalConsensus.subList(clippingBounds.getLeft(), clippingBounds.getRight() + 1);
            finalFastq = seqToFastq(fastq1.getReadName(), fastq1.getBaseQualityHeader(), clippedConsensus);
        }

        if(mDebugStatsWriter != null)
        {
            writeStats(fastq1, fastq2, hairpin1, hairpin2, rcMatch, seq1, seq2, seq1RC, seq2RC, forwardAlignment, forwardConsensus,
                    reverseAlignment, reverseConsensus, rcEndIndex, finalAlignment, finalConsensus, reverseConsensusTrimCount,
                    clippedConsensus, prefixTrimCount, suffixTrimCount);
        }

        return finalFastq;
    }

    private void writeStats(final FastqRecord fastq1, final FastqRecord fastq2, final HairpinInfo hairpin1, final HairpinInfo hairpin2,
            final ReverseComplementMatchInfo reverseComplementMatchInfo, final List<BaseQualPair> seq1, final List<BaseQualPair> seq2,
            final List<BaseQualPair> seq1RC, final List<BaseQualPair> seq2RC,
            @Nullable final List<Pair<BaseQualPair, BaseQualPair>> forwardAlignment,
            @Nullable final List<BaseQualPair> forwardConsensus, @Nullable final List<Pair<BaseQualPair, BaseQualPair>> reverseAlignment,
            @Nullable final List<BaseQualPair> reverseConsensus, int rcEndIndex,
            @Nullable final List<Pair<BaseQualPair, BaseQualPair>> finalAlignment,
            @Nullable final List<BaseQualPair> finalConsensus, int reverseConsensusTrimCount,
            @Nullable final List<BaseQualPair> clippedConsensus,
            int prefixTrimCount, int suffixTrimCount)
    {
        int trimmedLength = min(fastq1.getReadLength(), fastq2.getReadLength());
        int cutPoint = getCutPoint(trimmedLength, hairpin1, hairpin2, reverseComplementMatchInfo);

        // form initial output
        StringJoinerCounter statLine = new StringJoinerCounter(STAT_DELIMITER);
        statLine.add(fastq1.getReadName());
        statLine.add(String.valueOf(fastq1.getReadLength()));
        statLine.add(String.valueOf(fastq2.getReadLength()));
        statLine.add(fastq1.getReadString().substring(0, trimmedLength));
        statLine.add(sanatizeQualString(fastq1.getBaseQualityString().substring(0, trimmedLength)));
        statLine.add(fastq2.getReadString().substring(0, trimmedLength));
        statLine.add(sanatizeQualString(fastq2.getBaseQualityString().substring(0, trimmedLength)));

        HairpinInfo hairpin1_ = hairpin1 != null ? hairpin1 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin1_.StartIndex + 1));
        statLine.add(String.valueOf(hairpin1_.MatchCount));
        statLine.add(String.valueOf(hairpin1_.SuffixMatchLength));

        HairpinInfo hairpin2_ = hairpin2 != null ? hairpin2 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin2_.StartIndex + 1));
        statLine.add(String.valueOf(hairpin2_.MatchCount));
        statLine.add(String.valueOf(hairpin2_.SuffixMatchLength));

        ReverseComplementMatchInfo reverseComplementMatchInfo_ =
                reverseComplementMatchInfo != null ? reverseComplementMatchInfo : new ReverseComplementMatchInfo(0, -1, 0.0f, -1, 0.0f);
        statLine.add(String.valueOf(reverseComplementMatchInfo_.Read1Shift));
        statLine.add(String.valueOf(reverseComplementMatchInfo_.HighQualMismatchCount));
        statLine.add(format("%.4f", reverseComplementMatchInfo_.HighQualMismatchProportion));
        statLine.add(String.valueOf(reverseComplementMatchInfo_.TotalMismatchCount));
        statLine.add(format("%.4f", reverseComplementMatchInfo_.TotalMismatchProportion));

        if(cutPoint <= 0)
        {
            while(statLine.componentCount() < STAT_HEADERS.length)
            {
                statLine.add("-");
            }

            writeStatLine(mDebugStatsWriter, statLine.toString());
            return;
        }

        // forward naive consensus
        List<BaseQualPair> read1 = seq1.subList(0, cutPoint);
        List<BaseQualPair> read2 = seq2.subList(0, cutPoint);
        List<Pair<BaseQualPair, BaseQualPair>> naiveForwardAlignment = zip(read1, read2);
        FastqRecord naiveForwardConsensusFastq = seqToFastq(collapseAlignedSeqModC(naiveForwardAlignment));
        AlignmentStats naiveStats = getAlignmentStats(naiveForwardAlignment);

        // reverse naive consensus
        FastqRecord naiveReverseConsensusFastq = null;
        if(reverseComplementMatchInfo != null)
        {
            int rcStartIndex = -reverseComplementMatchInfo.Read1Shift;
            int rcLength = min(trimmedLength, cutPoint - rcStartIndex);
            if(rcLength > 0)
            {
                List<BaseQualPair> read1RC = seq1RC.subList(0, rcLength);
                List<BaseQualPair> read2RC = seq2RC.subList(0, rcLength);
                naiveReverseConsensusFastq = seqToFastq(collapseAlignedSeqModC(zip(read2RC, read1RC)));
            }
        }

        // forward aligned consensus
        AlignmentStats alignedStats = getAlignmentStats(forwardAlignment);
        FastqRecord forwardConsensusFastq = seqToFastq(forwardConsensus);

        int forwardMatchCount = 0;
        int forwardInsert1Count = 0;
        int forwardInsert2Count = 0;
        for(int i = 0; i < forwardAlignment.size(); i++)
        {
            BaseQualPair base1 = forwardAlignment.get(i).getLeft();
            BaseQualPair base2 = forwardAlignment.get(i).getRight();
            if(base1 != null && base2 != null)
            {
                forwardMatchCount++;
            }
            else if(base1 != null)
            {
                forwardInsert1Count++;
            }
            else
            {
                forwardInsert2Count++;
            }
        }

        // reverse aligned consensus
        int reverseMatchCount = 0;
        int reverseInsert1Count = 0;
        int reverseInsert2Count = 0;
        if(reverseAlignment != null)
        {
            for(int i = 0; i < reverseAlignment.size(); i++)
            {
                BaseQualPair base1 = reverseAlignment.get(i).getLeft();
                BaseQualPair base2 = reverseAlignment.get(i).getRight();
                if(base1 != null && base2 != null)
                {
                    reverseMatchCount++;
                }
                else if(base1 != null)
                {
                    reverseInsert1Count++;
                }
                else
                {
                    reverseInsert2Count++;
                }
            }
        }

        // final consensus
        int modCCMismatchCount = 0;
        int CmodCMismatchCount = 0;
        if(reverseAlignment != null)
        {
            for(Pair<BaseQualPair, BaseQualPair> alignedBase : finalAlignment)
            {
                BaseQualPair base1 = alignedBase.getLeft();
                BaseQualPair base2 = alignedBase.getRight();
                if(base1 == null || base2 == null)
                {
                    continue;
                }

                if(base1.Base == (byte) 'C' && base2.Base == MODC_BASE)
                {
                    CmodCMismatchCount++;
                }

                if(base2.Base == (byte) 'C' && base1.Base == MODC_BASE)
                {
                    modCCMismatchCount++;
                }
            }
        }

        // reverse resolved fastq comparison stats
        FastqRecord refFastq = mRefResolvesFastqMap.get(fastq1.getReadName());
        List<Pair<BaseQualPair, BaseQualPair>> refAlignment = null;
        int biomodalOffset = 0;
        int resolvedOverlappingBases = 0;
        int resolvedMissingMismatches = 0;
        int resolvedNonMissingMismatches = 0;
        int resolvedIndelCount = 0;
        if(refFastq != null && clippedConsensus != null)
        {
            List<BaseQualPair> refSeq = fastqToSeq(refFastq);
            refAlignment =
                    mAligner.align(clippedConsensus, refSeq, (a, b) -> getExactScore(a, b, true), OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, true, true, true, OptionalInt.empty());

            int firstMatchBase;
            for(firstMatchBase = 0; firstMatchBase < refAlignment.size(); firstMatchBase++)
            {
                Pair<BaseQualPair, BaseQualPair> alignedBase = refAlignment.get(firstMatchBase);
                BaseQualPair base1 = alignedBase.getLeft();
                BaseQualPair base2 = alignedBase.getRight();
                if(base1 != null && base2 != null)
                {
                    break;
                }

                biomodalOffset += base1 != null ? -1 : 1;
            }

            if(firstMatchBase == refAlignment.size())
            {
                refAlignment = null;
            }
            else
            {
                int lastMatchBase;
                for(lastMatchBase = refAlignment.size() - 1; lastMatchBase >= 0; lastMatchBase--)
                {
                    Pair<BaseQualPair, BaseQualPair> alignedBase = refAlignment.get(lastMatchBase);
                    BaseQualPair base1 = alignedBase.getLeft();
                    BaseQualPair base2 = alignedBase.getRight();
                    if(base1 != null && base2 != null)
                    {
                        break;
                    }
                }

                for(int i = firstMatchBase; i <= lastMatchBase; i++)
                {
                    Pair<BaseQualPair, BaseQualPair> alignedBase = refAlignment.get(i);
                    BaseQualPair base1 = alignedBase.getLeft();
                    BaseQualPair base2 = alignedBase.getRight();
                    if(base1 == null || base2 == null)
                    {
                        resolvedIndelCount++;
                        continue;
                    }

                    resolvedOverlappingBases++;

                    if(base1.Base != base2.Base)
                    {
                        if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE)
                        {
                            resolvedMissingMismatches++;
                        }
                        else
                        {
                            resolvedNonMissingMismatches++;
                        }
                    }
                }
            }
        }

        // form rest of output
        statLine.add(consensusReadForStatOutput(naiveForwardConsensusFastq.getReadString()));
        statLine.add(sanatizeQualString(naiveForwardConsensusFastq.getBaseQualityString()));

        statLine.add(String.valueOf(naiveStats.MissingCount));
        statLine.add(String.valueOf(naiveStats.MismatchCount));
        statLine.add(String.valueOf(naiveStats.HighQualMismatchCountGG));
        statLine.add(String.valueOf(naiveStats.HighQualMismatchCountOther));
        statLine.add(String.valueOf(naiveStats.LowQualUnambiguousCount));
        statLine.add(String.valueOf(naiveStats.LowQualAmbiguousCount));
        statLine.add(String.valueOf(naiveStats.ModCGCount));
        statLine.add(String.valueOf(naiveStats.ModCOtherCount));

        FastqRecord rcFastq1 = seqToFastq(seq1RC);
        FastqRecord rcFastq2 = seqToFastq(seq2RC);
        statLine.add(rcFastq1.getReadString());
        statLine.add(sanatizeQualString(rcFastq1.getBaseQualityString()));
        statLine.add(rcFastq2.getReadString());
        statLine.add(sanatizeQualString(rcFastq2.getBaseQualityString()));

        statLine.add(naiveReverseConsensusFastq == null ? "-" : consensusReadForStatOutput(naiveReverseConsensusFastq.getReadString()));
        statLine.add(naiveReverseConsensusFastq == null ? "-" : sanatizeQualString(naiveReverseConsensusFastq.getBaseQualityString()));

        FastqRecord forwardAlignmentRead1 = seqToFastq(getLeftElements(forwardAlignment, INS_BASE_QUAL_QUAL));
        FastqRecord forwardAlignmentRead2 = seqToFastq(getRightElements(forwardAlignment, INS_BASE_QUAL_QUAL));
        statLine.add(consensusReadForStatOutput(forwardAlignmentRead1.getReadString()));
        statLine.add(sanatizeQualString(forwardAlignmentRead1.getBaseQualityString()));
        statLine.add(consensusReadForStatOutput(forwardAlignmentRead2.getReadString()));
        statLine.add(sanatizeQualString(forwardAlignmentRead2.getBaseQualityString()));

        statLine.add(getCigar(forwardAlignment));
        statLine.add(String.valueOf(forwardMatchCount));
        statLine.add(String.valueOf(forwardInsert2Count));
        statLine.add(String.valueOf(forwardInsert1Count));
        statLine.add(String.valueOf(forwardInsert1Count + forwardInsert2Count));

        statLine.add(consensusReadForStatOutput(forwardConsensusFastq.getReadString()));
        statLine.add(sanatizeQualString(forwardConsensusFastq.getBaseQualityString()));

        statLine.add(String.valueOf(alignedStats.MissingCount));
        statLine.add(String.valueOf(alignedStats.MismatchCount));
        statLine.add(String.valueOf(alignedStats.HighQualMismatchCountGG));
        statLine.add(String.valueOf(alignedStats.HighQualMismatchCountOther));
        statLine.add(String.valueOf(alignedStats.LowQualUnambiguousCount));
        statLine.add(String.valueOf(alignedStats.LowQualAmbiguousCount));
        statLine.add(String.valueOf(alignedStats.ModCGCount));
        statLine.add(String.valueOf(alignedStats.ModCOtherCount));

        FastqRecord reverseAlignmentRead1 =
                reverseAlignment == null ? null : seqToFastq(getLeftElements(reverseAlignment, INS_BASE_QUAL_QUAL));
        FastqRecord reverseAlignmentRead2 =
                reverseAlignment == null ? null : seqToFastq(getRightElements(reverseAlignment, INS_BASE_QUAL_QUAL));
        statLine.add(reverseAlignment == null ? "-" : consensusReadForStatOutput(reverseAlignmentRead1.getReadString()));
        statLine.add(reverseAlignment == null ? "-" : sanatizeQualString(reverseAlignmentRead1.getBaseQualityString()));
        statLine.add(reverseAlignment == null ? "-" : consensusReadForStatOutput(reverseAlignmentRead2.getReadString()));
        statLine.add(reverseAlignment == null ? "-" : sanatizeQualString(reverseAlignmentRead2.getBaseQualityString()));
        statLine.add(reverseAlignment == null ? "-" : String.valueOf(rcEndIndex + 1));

        statLine.add(reverseAlignment == null ? "-" : getCigar(reverseAlignment));
        statLine.add(reverseAlignment == null ? "-" : String.valueOf(reverseMatchCount));
        statLine.add(reverseAlignment == null ? "-" : String.valueOf(reverseInsert2Count));
        statLine.add(reverseAlignment == null ? "-" : String.valueOf(reverseInsert1Count));
        statLine.add(reverseAlignment == null ? "-" : String.valueOf(reverseInsert1Count + reverseInsert2Count));

        FastqRecord reverseConsensusFastq = reverseAlignment == null ? null : seqToFastq(reverseConsensus);
        statLine.add(reverseAlignment == null ? "-" : consensusReadForStatOutput(reverseConsensusFastq.getReadString()));
        statLine.add(reverseAlignment == null ? "-" : sanatizeQualString(reverseConsensusFastq.getBaseQualityString()));

        FastqRecord finalAlignmentRead1 = seqToFastq(getLeftElements(finalAlignment, INS_BASE_QUAL_QUAL));
        FastqRecord finalAlignmentRead2 = seqToFastq(getRightElements(finalAlignment, INS_BASE_QUAL_QUAL));
        statLine.add(consensusReadForStatOutput(finalAlignmentRead1.getReadString()));
        statLine.add(sanatizeQualString(finalAlignmentRead1.getBaseQualityString()));
        statLine.add(consensusReadForStatOutput(finalAlignmentRead2.getReadString()));
        statLine.add(sanatizeQualString(finalAlignmentRead2.getBaseQualityString()));
        statLine.add(String.valueOf(reverseConsensusTrimCount));
        statLine.add(String.valueOf(modCCMismatchCount));
        statLine.add(String.valueOf(CmodCMismatchCount));

        FastqRecord finalConsensusFastq = seqToFastq(finalConsensus);
        statLine.add(consensusReadForStatOutput(finalConsensusFastq.getReadString()));
        statLine.add(sanatizeQualString(finalConsensusFastq.getBaseQualityString()));

        FastqRecord clippedConsensusFastq = clippedConsensus == null ? null : seqToFastq(clippedConsensus);
        statLine.add(String.valueOf(prefixTrimCount));
        statLine.add(String.valueOf(suffixTrimCount));
        statLine.add(clippedConsensus == null ? "-" : consensusReadForStatOutput(clippedConsensusFastq.getReadString()));
        statLine.add(clippedConsensus == null ? "-" : sanatizeQualString(clippedConsensusFastq.getBaseQualityString()));

        statLine.add(refFastq == null ? "-" : refFastq.getReadString());
        statLine.add(refFastq == null ? "-" : sanatizeQualString(refFastq.getBaseQualityString()));
        statLine.add(String.valueOf(clippedConsensus == null ? 0 : clippedConsensus.size()));
        statLine.add(String.valueOf(refFastq == null ? 0 : refFastq.getReadLength()));

        FastqRecord refAlignmentRead1 = refAlignment == null ? null : seqToFastq(getLeftElements(refAlignment, INS_BASE_QUAL_QUAL));
        FastqRecord refAlignmentRead2 = refAlignment == null ? null : seqToFastq(getRightElements(refAlignment, INS_BASE_QUAL_QUAL));
        statLine.add(refAlignment == null ? "-" : consensusReadForStatOutput(refAlignmentRead1.getReadString()));
        statLine.add(refAlignment == null ? "-" : sanatizeQualString(refAlignmentRead1.getBaseQualityString()));
        statLine.add(refAlignment == null ? "-" : consensusReadForStatOutput(refAlignmentRead2.getReadString()));
        statLine.add(refAlignment == null ? "-" : sanatizeQualString(refAlignmentRead2.getBaseQualityString()));
        statLine.add(refAlignment == null ? "-" : getCigar(refAlignment));

        statLine.add(String.valueOf(biomodalOffset));
        statLine.add(String.valueOf(resolvedOverlappingBases));
        statLine.add(String.valueOf(resolvedMissingMismatches));
        statLine.add(String.valueOf(resolvedNonMissingMismatches));
        statLine.add(String.valueOf(resolvedIndelCount));

        writeStatLine(mDebugStatsWriter, statLine.toString());
    }
}
