package com.hartwig.hmftools.linx.cn;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.cn.LohEvent.CN_DATA_NO_SV;
import static com.hartwig.hmftools.linx.types.SvVarData.NONE_SEGMENT_INFERRED;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.types.SvaConstants.SHORT_TI_LENGTH;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CnDataLoader
{
    private final String mPurpleDataPath;

    private int mRecordId;

    private DatabaseAccess mDbAccess;
    private List<StructuralVariantData> mSvDataList;
    private List<PurpleCopyNumber> mCnRecords;
    private Map<String,List<SvCNData>> mChrCnDataMap; // map of chromosome to CN data items
    private Map<Integer,SvCNData[]> mSvIdCnDataMap; // map of SV Ids to corresponding CN data pair
    private PurityContext mPurityContext;

    private List<LohEvent> mLohEventData;
    private List<HomLossEvent> mHomLossData;
    private Map<Integer,double[]> mSvPloidyCalcMap; // map of sample to SV Id & ploidy calc data
    private Map<String, double[]> mChrEndsCNMap; // telemore and centromere CN values

    public static double MIN_LOH_CN = 0.5;
    public static double TOTAL_CN_LOSS = 0.5;

    private static final Logger LOGGER = LogManager.getLogger(CnDataLoader.class);

    public CnDataLoader(final String purpleDataPath, DatabaseAccess dbAccess)
    {
        if(purpleDataPath.endsWith(File.separator))
            mPurpleDataPath = purpleDataPath;
        else
            mPurpleDataPath = purpleDataPath + File.separator;

        mDbAccess = dbAccess;
        mChrCnDataMap = Maps.newHashMap();
        mSvDataList = Lists.newArrayList();
        mCnRecords = Lists.newArrayList();
        mChrEndsCNMap = Maps.newHashMap();
        mSvIdCnDataMap = Maps.newHashMap();
        mSvPloidyCalcMap = Maps.newHashMap();
        mPurityContext = null;
        mLohEventData = Lists.newArrayList();
        mHomLossData = Lists.newArrayList();
        mRecordId = 0;
    }

    public final Map<Integer,double[]> getSvPloidyCalcMap() { return mSvPloidyCalcMap; }
    public final List<LohEvent> getLohData() { return mLohEventData; }
    public final List<HomLossEvent> getHomLossData() { return mHomLossData; }
    public final Map<String, double[]> getChrCopyNumberMap() { return mChrEndsCNMap; }
    public final Map<String,List<SvCNData>> getChrCnDataMap() { return mChrCnDataMap; }
    public final Map<Integer,SvCNData[]> getSvIdCnDataMap() { return mSvIdCnDataMap; }
    public final PurityContext getPurityContext() { return mPurityContext; }

    public void loadSampleData(final String sampleId, List<StructuralVariantData> svRecords)
    {
        mSvDataList.clear();
        mSvDataList.addAll(svRecords);

        loadCopyNumberData(sampleId);

        createChrCopyNumberMap();

        processSampleData(sampleId);
    }

    private void processSampleData(final String sampleId)
    {
        linkCopyNumberAndSvData(sampleId);

        findLohEvents(sampleId);

        reaclcAdjustedPloidy(sampleId);
    }

    public void clearState()
    {
        // shrink the data source to make future look-ups faster
        // and to release references to other SV objects (eg SvBreakend)
        mSvPloidyCalcMap.clear();
        mLohEventData.clear();
        mHomLossData.clear();
    }

    private void loadCopyNumberData(final String sampleId)
    {
        mChrCnDataMap.clear();
        mCnRecords.clear();

        if(mDbAccess == null)
        {
            try
            {
                mCnRecords = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(mPurpleDataPath, sampleId));
                mPurityContext = FittedPurityFile.read(mPurpleDataPath, sampleId);
            }
            catch(IOException e)
            {
                LOGGER.error("failed to load purity context: {}", e.toString());
                return;
            }

        }
        else
        {
            mCnRecords = mDbAccess.readCopynumbers(sampleId);
            LOGGER.debug("sample({}) retrievd {} CN entries", sampleId, mCnRecords.size());
            mPurityContext = mDbAccess.readPurityContext(sampleId);
        }

        String currentChromosome = "";
        List<SvCNData> cnDataList = null;
        for(final PurpleCopyNumber cnRecord : mCnRecords)
        {
            SvCNData cnData = new SvCNData(cnRecord, ++mRecordId);

            if(!currentChromosome.equals(cnData.Chromosome))
            {
                currentChromosome = cnData.Chromosome;
                cnDataList = Lists.newArrayList();
                mChrCnDataMap.put(cnData.Chromosome, cnDataList);
            }

            cnData.setIndex(cnDataList.size());
            cnDataList.add(cnData);
        }
    }

    private void linkCopyNumberAndSvData(final String sampleId)
    {
        mSvIdCnDataMap.clear();

        int unmatchedSVs = 0;

        for (final StructuralVariantData svData : mSvDataList)
        {
            if (svData.filter().equals(PON_FILTER_PON))
                continue;

            for(int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);

                final String svChromosome = isStart ? svData.startChromosome() : svData.endChromosome();
                long svPosition = isStart ? svData.startPosition() : svData.endPosition();
                long svOrientation = isStart ? svData.startOrientation() : svData.endOrientation();

                final List<SvCNData> cnDataList = mChrCnDataMap.get(svChromosome);

                if(cnDataList == null || cnDataList.isEmpty())
                    continue;

                boolean cnDataFound = false;

                for(SvCNData cnData : cnDataList)
                {
                    // ignore orientation during matching since CN change is not a reliable determinant of SV orientation
                    long cnPosition = cnData.StartPos;

                    if(svOrientation == 1)
                        --cnPosition;

                    if (svPosition != cnPosition)
                        continue;

                    if(!svData.type().toString().equals(cnData.SegStart))
                    {
                        if(svData.type() == SGL && svData.filter().equals(NONE_SEGMENT_INFERRED) && cnData.matchesSegment(NONE, true))
                        {
                            // SGL inferred == NONE in CN table
                        }
                        else
                        {
                            // allow matches if the position if correct
                            LOGGER.debug("SV({} chr={} pos={} type={}) matches {} CN segment: id({})",
                                    svData.id(), svChromosome, svPosition, svData.type(), cnData.SegStart, cnData.id());
                        }
                    }

                    // a match has been found
                    cnData.setStructuralVariantData(svData, isStart);

                    SvCNData[] cnDataItems = mSvIdCnDataMap.get(svData.id());

                    if(cnDataItems == null)
                    {
                        cnDataItems = new SvCNData[SE_END +1];
                        mSvIdCnDataMap.put(svData.id(), cnDataItems);
                    }

                    cnDataItems[be] = cnData;
                    cnDataFound = true;
                }

                if(!cnDataFound && svData.type() != INS)
                {
                    LOGGER.debug("SV({} chr={} pos={} type={}) {} unmatched)",
                            svData.id(), svChromosome, svPosition, svData.type(), isStart ? "start" : "end");

                    ++unmatchedSVs;
                }
            }
        }

        if(unmatchedSVs > 0)
        {
            LOGGER.warn("sample({}) has {} unmatched CN-SV segments", sampleId, unmatchedSVs);
        }
    }

    private StructuralVariantData getSvDataById(int svId)
    {
        for (final StructuralVariantData svData : mSvDataList)
        {
            if(svData.id() == svId)
                return svData;
        }

        return null;
    }

    private boolean isSingleVariant(final SvCNData cnData)
    {
        final StructuralVariantData svData = cnData.getStructuralVariantData();

        if(svData == null)
            return false;

        long svPosStart = svData.startOrientation() == -1 ? cnData.StartPos : cnData.StartPos - 1;
        long svPosEnd = cnData.EndPos + 1;

        return (svData.startPosition() == svPosStart && svData.endPosition() == svPosEnd);
    }

    private boolean skipMaleChromosomes(final String chromosome)
    {
        if(chromosome.equals("X") || chromosome.equals("Y"))
        {
            if (mPurityContext != null && mPurityContext.gender().toString().startsWith("MALE"))
            {
                return true;
            }
        }

        return false;
    }

    // private static String SPECIFIC_CHR = "9";
    private static String SPECIFIC_CHR = "";
    private static int REMOTE_SV_DISTANCE = 1000000;

    private void findLohEvents(final String sampleId)
    {
        mLohEventData.clear();
        mHomLossData.clear();

        // walk through the CN records looking for any loss of hetrozygosity, defined as a change in the actual baf to zero
        // and when CN rises back above this level, consider the section ended
        boolean isLohSection = false;
        double lohMinCN = 0;
        double lastMinCN = 0;
        double priorCN = 0; // before the LOH segment
        int lohSegments = 0;
        SvCNData lohStartCnData = null;
        int lohSectionCount = 0;
        int lohSVsMatchedCount = 0;
        boolean lohOnStartTelomere = false;
        boolean totalLoss = false;
        List<HomLossEvent> lohHomLossEvents = Lists.newArrayList();

        for(Map.Entry<String,List<SvCNData>> entry : mChrCnDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(skipMaleChromosomes(chromosome))
                continue;

            final List<SvCNData> cnDataList = entry.getValue();

            for(int index = 0; index < cnDataList.size(); ++index)
            {
                final SvCNData cnData = cnDataList.get(index);

                double minCN = cnData.minorAllelePloidy();

                boolean newChromosome = (index == 0);
                boolean reset = newChromosome;

                if (newChromosome && cnData.Chromosome.equals(SPECIFIC_CHR))
                {
                    LOGGER.debug("spec chr({})", SPECIFIC_CHR);
                }

                if (cnData.CopyNumber < TOTAL_CN_LOSS)
                {
                    // record homozygous loss
                    final SvCNData nextCnData = index < cnDataList.size() - 1 ? cnDataList.get(index + 1) : null;

                    StructuralVariantData svData = findSvData (cnData, 1);
                    StructuralVariantData nextSvData = nextCnData != null ? findSvData (nextCnData, -1) : null;

                    HomLossEvent homLoss = new HomLossEvent(chromosome, cnData.StartPos, cnData.EndPos,
                            cnData.SegStart, cnData.SegEnd,
                            svData != null ? svData.id() : CN_DATA_NO_SV,
                            nextSvData != null ? nextSvData.id() : CN_DATA_NO_SV);

                    LOGGER.debug("chromosome({}) HOM-loss at positions({} -> {}) segments({} - {})",
                            chromosome, homLoss.PosStart, homLoss.PosEnd, homLoss.SegStart, homLoss.SegEnd);

                    mHomLossData.add(homLoss);

                    if(isLohSection || lohOnStartTelomere)
                        lohHomLossEvents.add(homLoss);
                }

                if (isLohSection || lohOnStartTelomere)
                {
                    boolean lohRegained = (minCN >= MIN_LOH_CN);

                    // check that an SV with correct orientation exists here
                    if (lohRegained && cnData.matchesSV(true) && findSvData(cnData, -1) == null)
                    {
                        lohRegained = false;

                        LOGGER.debug(String.format("skipping LOH end segment(%s) with no SV match: type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                                cnData.asString(), cnData.SegStart,
                                cnData.ActualBaf, cnData.ObservedBaf, cnData.BafCount, cnData.CopyNumber));
                    }

                    if (lohRegained || reset)
                    {
                        // check for a short isolated TI and if found continue with the LOH
                        if (lohRegained && cnData.EndPos - cnData.StartPos <= SHORT_TI_LENGTH && index < cnDataList.size() - 1)
                        {
                            final SvCNData nextData = cnDataList.get(index + 1);
                            double nextMinCN = nextData.minorAllelePloidy();

                            if (nextData.EndPos - cnData.StartPos > REMOTE_SV_DISTANCE && nextMinCN < MIN_LOH_CN
                            && lohStartCnData != null && cnData.StartPos - lohStartCnData.StartPos > REMOTE_SV_DISTANCE)
                            {
                                LOGGER.debug("chr({}) skipping short isolated TI(id={} {} pos={} length={})",
                                        chromosome, cnData.id(), cnData.SegStart, cnData.StartPos, cnData.EndPos - cnData.StartPos);
                                continue;
                            }
                        }

                        if (lohOnStartTelomere || totalLoss)
                        {
                            // LOH section invalidated
                            processLOHData(chromosome, lohStartCnData, cnData, priorCN, lohMinCN,
                                    lohSegments, false, lohHomLossEvents);

                            lohOnStartTelomere = false;
                            totalLoss = false;
                        }
                        else
                        {
                            // log all relevant data for this completed section
                            lohSVsMatchedCount += processLOHData(chromosome, lohStartCnData, cnData, priorCN, lohMinCN,
                                    lohSegments, false, lohHomLossEvents);
                            ++lohSectionCount;
                        }

                        reset = true;
                    }
                    else if (cnData.matchesSegment(TELOMERE, false))
                    {
                        // rest of arm was lost so no linking SV for LOH section - but still record the event
                        processLOHData(chromosome, lohStartCnData, cnData, priorCN, lohMinCN,
                                lohSegments, true, lohHomLossEvents);
                        reset = true;
                    }
                    else if (cnData.CopyNumber < TOTAL_CN_LOSS)
                    {
                        // other chromatid loss has occurred - this will cancel a valid LOH due to uncertainty unless due to a simple SV
                        if (cnData.matchesSegment(SegmentSupport.DEL, true) && cnData.matchesSegment(SegmentSupport.DEL, false)
                                && isSingleVariant(cnData))
                        {
                            LOGGER.debug("total CN loss matches single SV({} : {} -> {})", chromosome, cnData.StartPos, cnData.EndPos);
                        }
                        else
                        {
                            totalLoss = true;
                        }
                    }
                    else
                    {
                        ++lohSegments;
                        lohMinCN = min(lohMinCN, minCN);
                    }
                }
                else if (!reset)
                {
                    boolean lohLost = minCN < MIN_LOH_CN;

                    // check that an SV with correct orientation exists here
                    if (lohLost && cnData.matchesSV(true) && findSvData(cnData, 1) == null)
                    {
                        lohLost = false;
                    }

                    if (lohLost)
                    {
                        // new LOH section identified
                        isLohSection = true;
                        lohSegments = 1;
                        lohStartCnData = cnData;
                        lohMinCN = minCN;
                        priorCN = lastMinCN;

                        if (lohOnStartTelomere)
                        {
                            LOGGER.debug("chr({}) LOH at telemore", chromosome);
                        }
                        else
                        {
                            LOGGER.debug(String.format("chr(%s) starting LOH at pos(%d) minCN(%.3f cn=%.3f baf=%.3f) priorCN(%.3f)",
                                    chromosome, cnData.StartPos, cnData.CopyNumber, cnData.ActualBaf, lohMinCN, priorCN));

                            // check for segments ending on telomere
                            if (cnData.matchesSegment(TELOMERE, false))
                            {
                                processLOHData(chromosome, lohStartCnData, lohStartCnData, priorCN, lohMinCN,
                                        lohSegments, true, null);
                                reset = true;
                            }
                        }
                    }
                }

                if (reset)
                {
                    isLohSection = false;
                    lohOnStartTelomere = false;
                    lohSegments = 0;
                    lohStartCnData = null;
                    totalLoss = false;
                    lohHomLossEvents.clear();

                    if (newChromosome && (minCN < MIN_LOH_CN))
                    {
                        lohOnStartTelomere = true;
                        lohStartCnData = cnData;
                    }
                }

                lastMinCN = minCN;
            }
        }

        LOGGER.debug("sample({}) LOH sections({}) fullMatched({})",
                sampleId, lohSectionCount, lohSVsMatchedCount);
    }

    private StructuralVariantData findSvData(final SvCNData cnData, int requiredOrient)
    {
        if (cnData.matchesSegment(UNKNOWN, true)
        || cnData.matchesSegment(TELOMERE, true)
        || cnData.matchesSegment(CENTROMERE, true))
        {
            return null;
        }

        long svPosition = requiredOrient == -1 ? cnData.StartPos : cnData.StartPos - 1;

        final StructuralVariantData svData = cnData.getStructuralVariantData();

        if (svData == null)
            return null;

        if (cnData.svLinkOnStart() && svData.startOrientation() == requiredOrient && svData.startPosition() == svPosition)
            return svData;
        else if (!cnData.svLinkOnStart() && svData.endOrientation() == requiredOrient && svData.endPosition() == svPosition)
            return svData;
        else
            return null;
    }

    private int processLOHData(final String chr, SvCNData startData, SvCNData endData,
            double lastMinCN, double lohMinCN, int segCount, boolean incomplete, List<HomLossEvent> lohHomLossEvents)
    {
        StructuralVariantData startSvData = findSvData(startData, 1);

        StructuralVariantData endSvData = null;
        long lohLength = 0;

        if(incomplete)
        {
            // ended at the telomere
            lohLength = endData.EndPos - startData.StartPos;
            endSvData = null;
        }
        else if(endData.Chromosome.equals(chr) && !startData.matchesSegment(TELOMERE, false))
        {
            lohLength = endData.StartPos - startData.StartPos;
            endSvData = findSvData(endData, -1);
        }
        else
        {
            // segment has either started and/or finished on the telomere segment
            endData = startData;
            lohLength = startData.EndPos - startData.StartPos + 1;
        }

        if(LOGGER.isDebugEnabled())
        {
            if (lohLength <= 0)
            {
                LOGGER.warn("negative length({})", lohLength);
            }

            if (startSvData != null && endSvData != null)
            {
                if (startSvData.id() == endSvData.id())
                {
                    LOGGER.debug("cnID({} -> {}) matches singleSV({} - {})",
                            startData.asString(), endData.asString(), startSvData.id(), startSvData.type());
                }
                else
                {
                    LOGGER.debug("cnID({} -> {}) matches pairSV({} -> {})",
                            startData.asString(), endData.asString(), startSvData.id(), endSvData.id());
                }
            }
            else
            {
                if (!incomplete && lohHomLossEvents == null)
                {
                    if (startSvData == null && startData.matchesSV(true))
                    {
                        boolean incorrectOrientation = (startData.getStructuralVariantData() != null);

                        LOGGER.debug(String.format("LOH start segment(%s) no SV match: orient(%s) type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                                startData.asString(), incorrectOrientation ? "wrong" : "ok", startData.SegStart,
                                startData.ActualBaf, startData.ObservedBaf, startData.BafCount, startData.CopyNumber));
                    }

                    if (!incomplete && endSvData == null && endData.matchesSV(true))
                    {
                        boolean incorrectOrientation = (endData.getStructuralVariantData() != null);

                        LOGGER.debug(String.format("LOH end segment(%s) no SV match: orient(%s), type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                                endData.asString(), incorrectOrientation ? "wrong" : "ok", endData.SegStart,
                                endData.ActualBaf, endData.ObservedBaf, endData.BafCount, endData.CopyNumber));
                    }
                }

                LOGGER.debug("cnID({} & {}) not fully matched pairSV({} -> {})",
                        startData.asString(), startData.SegStart, endData.asString(), endData.SegStart,
                        startSvData != null ? startSvData.id() : "", endSvData != null ? endSvData.id() : "");
            }
        }

        LohEvent lohData = new LohEvent(
                chr,
                startData.StartPos, incomplete ? endData.EndPos : endData.StartPos,
                startData.SegStart, incomplete ? endData.SegEnd : endData.SegStart,
                lastMinCN, startData.minorAllelePloidy(), endData.minorAllelePloidy(),
                lohMinCN, segCount, lohLength,
                startSvData != null ? startSvData.id() : CN_DATA_NO_SV,
                endSvData != null ? endSvData.id() : CN_DATA_NO_SV);

        if(lohHomLossEvents != null)
            lohData.addHomLossEvents(lohHomLossEvents);

        mLohEventData.add(lohData);

        return (startSvData != null && endSvData != null) ? 1 : 0;
    }


    public static int P_ARM_TELOMERE_CN = 0;
    public static int CENTROMERE_CN = 1;
    public static int Q_ARM_TELOMERE_CN = 2;

    public void createChrCopyNumberMap()
    {
        mChrEndsCNMap.clear();

        double[] currentChrData = null;

        for(final PurpleCopyNumber cnRecord : mCnRecords)
        {
            if(cnRecord.segmentStartSupport().equals(TELOMERE))
            {
                currentChrData = new double[Q_ARM_TELOMERE_CN+1];
                mChrEndsCNMap.put(cnRecord.chromosome(), currentChrData);
                currentChrData[P_ARM_TELOMERE_CN] = cnRecord.averageTumorCopyNumber();
            }
            else if(cnRecord.segmentStartSupport().equals(CENTROMERE))
            {
                currentChrData[CENTROMERE_CN] = cnRecord.averageTumorCopyNumber();
            }
            else if(cnRecord.segmentEndSupport().equals(TELOMERE))
            {
                currentChrData[Q_ARM_TELOMERE_CN] = cnRecord.averageTumorCopyNumber();
            }
        }
    }

    public static int CN_SEG_DATA_CN_BEFORE = 0;
    public static int CN_SEG_DATA_MAP_BEFORE = 1;
    public static int CN_SEG_DATA_CN_AFTER = 2;
    public static int CN_SEG_DATA_MAP_AFTER = 3;

    public double[] getCentromereCopyNumberData(final String chromosome, boolean fromPArm)
    {
        // get copy number and major allele ploidy for the segment leading up and leading from the centromere
        List<SvCNData> cnDataList = mChrCnDataMap.get(chromosome);

        double[] results = new double[CN_SEG_DATA_MAP_AFTER+1];

        if(cnDataList == null || cnDataList.isEmpty())
            return results;

        for(int i = 0; i < cnDataList.size(); ++i)
        {
            SvCNData cnData = cnDataList.get(i);

            if(cnData.matchesSegment(CENTROMERE, false))
            {
                SvCNData nextCnData = cnDataList.get(i+1);

                if(fromPArm)
                {
                    results[CN_SEG_DATA_CN_BEFORE] = cnData.CopyNumber;
                    results[CN_SEG_DATA_MAP_BEFORE] = cnData.majorAllelePloidy();
                    results[CN_SEG_DATA_CN_AFTER] = nextCnData.CopyNumber;
                    results[CN_SEG_DATA_MAP_AFTER] = nextCnData.majorAllelePloidy();
                }
                else
                {
                    results[CN_SEG_DATA_CN_BEFORE] = nextCnData.CopyNumber;
                    results[CN_SEG_DATA_MAP_BEFORE] = nextCnData.majorAllelePloidy();
                    results[CN_SEG_DATA_CN_AFTER] = cnData.CopyNumber;
                    results[CN_SEG_DATA_MAP_AFTER] = cnData.majorAllelePloidy();
                }

                break;
            }
        }

        return results;
    }


    private void reaclcAdjustedPloidy(final String sampleId)
    {
        mSvPloidyCalcMap.clear();

        for (Map.Entry<Integer, SvCNData[]> entry : mSvIdCnDataMap.entrySet())
        {
            final int svId = entry.getKey();
            final SvCNData[] cnDataPair = entry.getValue();

            final SvCNData cnStartData = cnDataPair[SE_START];

            if (cnStartData == null)
            {
                LOGGER.error("SV({}) missing start copy number data", svId);
                continue;
            }

            StructuralVariantData svData = cnStartData.getStructuralVariantData();

            if(svData.id() != svId)
            {
                svData = getSvDataById(svId);
            }

            final SvCNData cnStartPrevData = getCNSegment(cnStartData.Chromosome,  cnStartData.getIndex() - 1);

            final SvCNData cnEndNextData = cnDataPair[SE_END];
            final SvCNData cnEndData = cnEndNextData != null ? getCNSegment(cnEndNextData.Chromosome, cnEndNextData.getIndex() - 1) : null;

            int tumorReadCountStart = svData.startTumorVariantFragmentCount();
            int tumorReadCountEnd = svData.endTumorVariantFragmentCount();

            double adjVafStart = svData.adjustedStartAF();
            double adjVafEnd = svData.adjustedEndAF();

            double maxCNStart = max(cnStartPrevData.CopyNumber, cnStartData.CopyNumber);
            double maxCNEnd = cnEndData != null && cnEndNextData != null ? max(cnEndData.CopyNumber, cnEndNextData.CopyNumber): 0;

            final int[] startDepthData = { cnStartPrevData.DepthWindowCount, cnStartData.DepthWindowCount };
            double cnChgStart = svData.adjustedStartCopyNumberChange();

            final int[] endDepthData = { cnEndData != null ? cnEndData.DepthWindowCount : 0,
                    cnEndNextData != null ? cnEndNextData.DepthWindowCount : 0 };

            double cnChgEnd = svData.adjustedEndCopyNumberChange();

            // for the special case of a foldback inversion, the CN change isn't reliable for the breakends, but the CN change over the region is reliable
            // so use this value instead
            if(cnEndData != null && svData.type() == INV && (cnStartData.getIndex() == cnEndData.getIndex()))
            {
                if(cnStartPrevData.DepthWindowCount > cnStartData.DepthWindowCount && cnEndNextData.DepthWindowCount > cnStartData.DepthWindowCount)
                {
                    double cnChange = abs(cnStartPrevData.CopyNumber - cnEndNextData.CopyNumber) * 0.5;

                    // use the DWC from the outer segments
                    endDepthData[0] = endDepthData[1];
                    startDepthData[1] = startDepthData[0];
                    cnChgStart = cnChgEnd = cnChange;
                }
            }

            final double calcResults[] = calcAdjustedPloidyValues(cnChgStart, cnChgEnd,
                    tumorReadCountStart, tumorReadCountEnd, adjVafStart, adjVafEnd, maxCNStart, maxCNEnd,
                    startDepthData, cnEndData != null ? endDepthData : null);

            double ploidyEstimate = calcResults[APC_EST_PLOIDY];
            double ploidyUncertainty = calcResults[APC_EST_UNCERTAINTY];

            if(ploidyUncertainty == 0 || Double.isNaN(ploidyEstimate) || Double.isNaN(ploidyUncertainty))
            {
                LOGGER.debug("sample({}) svID({} type={}) unexpected ploidy(est={} unc={})",
                        sampleId, svData.id(), svData.type(), ploidyEstimate, ploidyUncertainty);
            }

            if(Double.isNaN(ploidyEstimate))
                ploidyEstimate = 0;

            if(Double.isNaN(ploidyUncertainty))
                ploidyUncertainty = 0;

            double[] ploidyCalcs = {ploidyEstimate, ploidyUncertainty};
            mSvPloidyCalcMap.put(svData.id(), ploidyCalcs);

        }
    }

    public static int APC_EST_PLOIDY = 0;
    public static int APC_EST_UNCERTAINTY = 1;
    public static int APC_RESULTS_MAX = 2;

    private static double POIS_PROB_LOW = 0.005;
    private static double POIS_PROB_HIGH = 0.995;

    public static double[] calcAdjustedPloidyValues(double cnChgStart, double cnChgEnd,
            int tumorReadCountStart, int tumorReadCountEnd, double adjVafStart, double adjVafEnd,
            double maxCNStart, double maxCNEnd, final int[] startDepthData, final int[] endDepthData)
    {
        double cnUncertaintyStart = calcCopyNumberSideUncertainty(maxCNStart, startDepthData);
        double cnUncertaintyEnd = endDepthData != null ? calcCopyNumberSideUncertainty(maxCNEnd, endDepthData) : 0;

        List<Double> observations = Lists.newArrayList();
        List<Double> uncertainties = Lists.newArrayList();

        if(cnUncertaintyStart > 0)
        {
            observations.add(cnChgStart);
            uncertainties.add(cnUncertaintyStart);
        }

        if(cnUncertaintyEnd > 0)
        {
            observations.add(cnChgEnd);
            uncertainties.add(cnUncertaintyEnd);
        }

        // calculate a ploidy for each side independently
        if(adjVafStart > 0)
        {
            double ploidyStart = adjVafStart * maxCNStart;
            double ploidyUncertaintyStart = calcPloidyUncertainty(tumorReadCountStart, adjVafStart, maxCNStart, cnUncertaintyStart);

            if (ploidyUncertaintyStart > 0)
            {
                observations.add(ploidyStart);
                uncertainties.add(ploidyUncertaintyStart);
            }
        }

        if(adjVafEnd > 0)
        {
            double ploidyEnd = adjVafEnd * maxCNEnd;
            double ploidyUncertaintyEnd = calcPloidyUncertainty(tumorReadCountEnd, adjVafEnd, maxCNEnd, cnUncertaintyEnd);

            if (ploidyUncertaintyEnd > 0)
            {
                observations.add(ploidyEnd);
                uncertainties.add(ploidyUncertaintyEnd);
            }
        }

        double estPloidy = 0;
        double estUncertainty = 0;

        if(observations.size() == 1)
        {
            estPloidy = observations.get(0);
            estUncertainty = uncertainties.get(0);
        }
        else
        {
            double sumUncertainty = 0;
            double sumObservedUncertainty = 0;

            for(int i = 0; i < observations.size(); ++i)
            {
                double uncertInvSqrd = 1 / pow(uncertainties.get(i), 2);
                sumUncertainty += uncertInvSqrd;
                sumObservedUncertainty += observations.get(i) * uncertInvSqrd;
            }

            // consloidatedPloidy =  SUM[Observation(i)*(1/Uncertainty(i)^2)] / Sum[1/Uncertainty(i)^2]
            estPloidy = sumObservedUncertainty / sumUncertainty;

            double adjUncertainty = 0;

            for(int i = 0; i < observations.size(); ++i)
            {
                double uncertInvSqrd = 1 / pow(uncertainties.get(i), 2);
                double relativeUncertainty = uncertInvSqrd * pow(max(observations.get(i) - estPloidy, uncertainties.get(i)/2),2);
                adjUncertainty += relativeUncertainty;
            }

            // consolidatedUncertainty = SQRT(countObservations/(countObervations-1) * SUM[(1/Uncertainty(i)^2*(MAX(Observation(i)-consolidatedPloidy,Uncertainty(i)/2))^2]
            // / Sum[1/Uncertainty(i)^2] )

            estUncertainty = sqrt(observations.size() / (observations.size() - 1) * adjUncertainty / sumUncertainty);
        }

        // populate a results object
        double[] results = new double[APC_RESULTS_MAX+1];

        results[APC_EST_UNCERTAINTY] = estUncertainty;
        results[APC_EST_PLOIDY] = estPloidy;

        return results;
    }

    private static double ABS_UNCERTAINTY = 0.15;
    private static double RELATIVE_UNCERTAINTY = 0.1;
    private static double ADDITIONAL_ABS_UNCERTAINTY = 0.4;
    private static double ADDITIONAL_REL_UNCERTAINTY = 0.15;
    private static double PROPORTION_CNCHANGE_USED_IN_PLOIDY_UNC = 0.5;
    private static double NO_DEPTH_CNCHANGE_UNC = 0.5;

    private static double calcPloidyUncertainty(int tumorReadCount, double adjVaf, double maxCopyNumber, double cnChangeUncertainty)
    {
        if(cnChangeUncertainty <= 0)
            return 0;

        double poissonRCLow = calcPoisonReadCount(tumorReadCount, POIS_PROB_LOW);
        double poissonRCHigh = calcPoisonReadCount(tumorReadCount, POIS_PROB_HIGH);

        double poissonVafLow = adjVaf * poissonRCLow / tumorReadCount;
        double poissonVafHigh = adjVaf * poissonRCHigh / tumorReadCount;

        double ploidyUncertainty = (poissonVafHigh * (maxCopyNumber + PROPORTION_CNCHANGE_USED_IN_PLOIDY_UNC * cnChangeUncertainty));
        ploidyUncertainty -= (poissonVafLow * (maxCopyNumber - PROPORTION_CNCHANGE_USED_IN_PLOIDY_UNC * cnChangeUncertainty));
        ploidyUncertainty *= 0.5;

        return ploidyUncertainty;
    }

    private static double calcCopyNumberSideUncertainty(double copyNumber, final int[] depthData)
    {
        int minDepthCount = min(depthData[0], depthData[1]);

        if(minDepthCount <= 0)
            return NO_DEPTH_CNCHANGE_UNC * copyNumber;

        double uncertainty = max(copyNumber*RELATIVE_UNCERTAINTY, ABS_UNCERTAINTY);
        uncertainty += max(ADDITIONAL_ABS_UNCERTAINTY, ADDITIONAL_REL_UNCERTAINTY * copyNumber) / sqrt(minDepthCount);

        return uncertainty;
    }

    private static double calcPoisonReadCount(int readCount, double requiredProb)
    {
        return calcPoissonObservedGivenProb(readCount, requiredProb);
    }

    private static int calcPoissonObservedGivenProb(int expectedVal, double requiredProb)
    {
        if(expectedVal <= 0)
            return 0;

        PoissonDistribution poisson = new PoissonDistribution(expectedVal);

        int maxIterations = 20;
        int iterations = 0;

        double refCount = 25;
        double refRangePerc = 0.44;
        double rangePerc = refRangePerc / sqrt(expectedVal/refCount);
        double range = rangePerc * expectedVal;

        int testValueUpper;
        int testValueLower;

        if(requiredProb > 0.5)
        {
            testValueUpper = (int) round(expectedVal + range * 1.5);
            testValueLower = (int) round(max(expectedVal + range * 0.2, 0));
        }
        else
        {
            testValueUpper = (int) round(expectedVal - range * 0.2);
            testValueLower = (int) round(max(expectedVal - range * 1.2, 0));
        }

        int testValue = (int)((testValueLower + testValueUpper) * 0.5);

        double currentProb = poisson.cumulativeProbability(testValue);
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff < 0.001)
                break;

            // if prob is too high, need to lower the test value
            if(currentProb > requiredProb)
            {
                if(testValue <= testValueLower + 1)
                    break;

                testValueUpper = testValue;
                testValue = (int)round((testValue + testValueLower) * 0.5);
            }
            else
            {
                if(testValue >= testValueUpper - 1)
                    break;

                testValueLower = testValue;
                testValue = (int)round((testValue + testValueUpper) * 0.5);
            }

            currentProb = poisson.cumulativeProbability(testValue);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            LOGGER.warn(String.format("max iterations reached: value(%d) test(%d) prob(%.4f diff=%.4f)",
                    expectedVal, testValue, currentProb, probDiff));
        }

        return testValue;
    }

    public final SvCNData getCNSegment(final String chromosome, int index)
    {
        final List<SvCNData> cnDataList = mChrCnDataMap.get(chromosome);

        if(cnDataList == null || cnDataList.isEmpty())
            return null;

        if(index < 0 || index >= cnDataList.size())
            return null;

        return cnDataList.get(index);
    }

}
