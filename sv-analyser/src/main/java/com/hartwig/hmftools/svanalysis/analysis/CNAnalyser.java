package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.NONE_SEGMENT_INFERRED;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvLOH;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CNAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(CNAnalyser.class);

    private final String mOutputPath;

    private List<String> mSampleIds;
    private int mRecordId;
    private int mNoneSvId;

    private boolean mWriteAdjustedPloidyToFile;
    private boolean mUpdateAdjustedPloidyToDB;
    private boolean mWriteVerbosePloidyData;

    BufferedWriter mFileWriter;
    BufferedWriter mRecalcPloidyFileWriter;
    DatabaseAccess mDbAccess;
    List<StructuralVariantData> mSvDataList;
    List<SvCNData> mCnDataList;
    Map<String,SvCNData[]> mSvCnDataMap;

    Map<String, List<SvLOH>> mSampleLohData;
    private Map<String,Map<String,double[]>> mSampleSvPloidyCalcMap;

    public static final String SV_PLOIDY_CALC_FILE = "sv_ploidy_file";
    private static final String WRITE_PLOIDY_TO_FILE = "write_ploidy_to_file";
    private static final String UPDATE_PLOIDY_TO_DB = "update_ploidy_to_db";
    private static final String WRITE_VERBOSE_PLOIDY_DATA = "verbose_ploidy_data";

    public static double MIN_LOH_CN = 0.5;
    public static double TOTAL_CN_LOSS = 0.5;

    public CNAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mOutputPath = outputPath;
        mSampleIds = Lists.newArrayList();
        mFileWriter = null;
        mRecalcPloidyFileWriter = null;
        mWriteAdjustedPloidyToFile = false;
        mWriteVerbosePloidyData = false;
        mDbAccess = dbAccess;
        mCnDataList = Lists.newArrayList();
        mSvDataList = Lists.newArrayList();
        mSvCnDataMap = new HashMap();
        mSampleSvPloidyCalcMap = new HashMap();
        mSampleLohData = null;
        mUpdateAdjustedPloidyToDB = false;
        mRecordId = 0;
        mNoneSvId = 0;
    }

    public final Map<String,Map<String,double[]>> getSampleSvPloidyCalcMap() { return mSampleSvPloidyCalcMap; }
    public final Map<String, List<SvLOH>> getSampleLohData() { return mSampleLohData; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SV_PLOIDY_CALC_FILE, true, "SV_PLOIDY_CALC_FILE");
        options.addOption(WRITE_PLOIDY_TO_FILE, false, "Write adjusted ploidy to CSV");
        options.addOption(WRITE_VERBOSE_PLOIDY_DATA, false, "Write all ploidy calc working data");
        options.addOption(UPDATE_PLOIDY_TO_DB, false, "Update SV table with ploidy calcs");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> sampleIds)
    {
        mSampleIds.addAll(sampleIds);
        mWriteAdjustedPloidyToFile = cmd.hasOption(WRITE_PLOIDY_TO_FILE);
        mWriteVerbosePloidyData = cmd.hasOption(WRITE_VERBOSE_PLOIDY_DATA);
        mUpdateAdjustedPloidyToDB = cmd.hasOption(UPDATE_PLOIDY_TO_DB);
        return true;
    }

    public void run()
    {
        if(mDbAccess == null)
            return;

        int sampleCount = 0;

        for(final String sampleId : mSampleIds)
        {
            loadCopyNumberData(sampleId);

            ++sampleCount;
            LOGGER.info("analysing sample({}) with {} CN entries, totalProcessed({})", sampleId, mCnDataList.size(), sampleCount);

            loadSVData(sampleId);
            linkCopyNumberAndSvData(sampleId);

            findLohEvents(sampleId);

            if(mWriteAdjustedPloidyToFile || mUpdateAdjustedPloidyToDB)
                reaclcAdjustedPloidy(sampleId);
        }
    }

    private void loadCopyNumberData(final String sampleId)
    {
        mCnDataList.clear();
        List<PurpleCopyNumber> cnRecords = mDbAccess.readCopynumbers(sampleId);

        for(final PurpleCopyNumber cnRecord : cnRecords)
        {
            SvCNData cnData = new SvCNData(cnRecord, ++mRecordId);
            cnData.setIndex(mCnDataList.size());
            mCnDataList.add(cnData);
        }
    }

    private void loadSVData(final String sampleId)
    {
        if(mDbAccess == null)
            return;

        mSvDataList = mDbAccess.readStructuralVariantData(sampleId);

        if(mSvDataList.isEmpty())
        {
            LOGGER.warn("sample({}) no SV records found", sampleId);
        }
    }

    private void linkCopyNumberAndSvData(final String sampleId)
    {
        mSvCnDataMap.clear();

        int unmatchedSVs = 0;
        for(SvCNData cnData : mCnDataList)
        {
            if (!cnData.matchesSV(true) && !cnData.matchesSegment(NONE, true))
                continue;

            // ignore orientation during matching since CN change is not a reliable determinant of SV orientation
            long cnPosition = cnData.StartPos;

            for (final StructuralVariantData svData : mSvDataList)
            {
                if (svData.filter().equals(PON_FILTER_PON))
                    continue;

                if(!svData.type().toString().equals(cnData.SegStart))
                {
                    if(svData.type() == SGL && svData.filter().equals(NONE_SEGMENT_INFERRED) && cnData.matchesSegment(NONE, true))
                    {
                        // SGL inferred == NONE in CN table
                    }
                    else if(cnData.matchesSegment(MULTIPLE, true))
                    {
                        // also valid
                    }
                    else
                    {
                        continue;
                    }
                }

                long svPosStart = svData.startOrientation() == -1 ? cnPosition : cnPosition - 1;
                long svPosEnd = svData.endOrientation() == -1 ? cnPosition : cnPosition - 1;
                boolean matchOnStart = false;

                if (svData.startChromosome().equals(cnData.Chromosome) && svData.startPosition() == svPosStart)
                {
                    matchOnStart = true;
                }
                else if (svData.endChromosome().equals(cnData.Chromosome) && svData.endPosition() == svPosEnd)
                {
                    matchOnStart = false;
                }
                else
                {
                    continue;
                }

                cnData.setStructuralVariantData(svData, matchOnStart);

                SvCNData[] cnDataItems = mSvCnDataMap.get(svData.id());

                if(cnDataItems == null)
                {
                    cnDataItems = new SvCNData[2];
                    mSvCnDataMap.put(svData.id(), cnDataItems);
                }

                cnDataItems[matchOnStart ? SVI_START : SVI_END] = cnData;
            }
        }

        if(unmatchedSVs > 0)
        {
            LOGGER.warn("sample({}) has {} unmatched CN-SV segments", sampleId, unmatchedSVs);
        }
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

    // private static String SPECIFIC_CHR = "12";
    private static String SPECIFIC_CHR = "";
    private static int REMOTE_SV_DISTANCE = 1000000;

    private void findLohEvents(final String sampleId)
    {
        String currentChr = "";

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
        boolean totalLoss = false;StructuralVariantData totalLossSv = null;

        for(int index = 0; index < mCnDataList.size(); ++index)
        {
            final SvCNData cnData = mCnDataList.get(index);

            double minCN = (1 - cnData.ActualBaf) * cnData.CopyNumber;

            boolean newChromosome = currentChr.isEmpty() || (!currentChr.isEmpty() && !cnData.Chromosome.equals(currentChr));
            boolean reset = newChromosome;

            if(newChromosome && cnData.Chromosome.equals(SPECIFIC_CHR))
            {
                LOGGER.debug("spec chr({})", SPECIFIC_CHR);
            }

            if(isLohSection || lohOnStartTelomere)
            {
                boolean lohRegained = (minCN >= MIN_LOH_CN);

                // check that an SV with correct orientation exists here
                if(lohRegained && cnData.matchesSV(true) && findSvData(cnData, -1) == null)
                {
                    lohRegained = false;

                    LOGGER.debug(String.format("skipping LOH end segment(%s) with no SV match: type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                            cnData.asString(), cnData.SegStart,
                            cnData.ActualBaf, cnData.ObservedBaf, cnData.BafCount, cnData.CopyNumber));
                }

                if(lohRegained || reset)
                {
                    // check for a short isolated TI and if found continue with the LOH
                    if(lohRegained && cnData.EndPos - cnData.StartPos <= SHORT_TI_LENGTH && index < mCnDataList.size() - 1)
                    {
                        final SvCNData nextData = mCnDataList.get(index+1);
                        double nextMinCN = (1 - nextData.ActualBaf) * nextData.CopyNumber;

                        if(nextData.EndPos - cnData.StartPos > REMOTE_SV_DISTANCE && nextMinCN < MIN_LOH_CN
                        && lohStartCnData!= null && cnData.StartPos - lohStartCnData.StartPos > REMOTE_SV_DISTANCE)
                        {
                            LOGGER.debug("chr({}) skipping short isolated TI(id={} {} pos={} length={})",
                                    currentChr, cnData.id(), cnData.SegStart, cnData.StartPos, cnData.EndPos - cnData.StartPos);

                            writeLOHData(sampleId, currentChr, cnData, nextData, priorCN, lohMinCN,
                                    1, false, true, false);
                            continue;
                        }
                    }

                    if(lohOnStartTelomere || totalLoss)
                    {
                        // LOH section invalidated
                        writeLOHData(sampleId, currentChr, lohStartCnData, cnData, priorCN, lohMinCN,
                                lohSegments, false, false, !totalLoss);

                        lohOnStartTelomere = false;
                        totalLoss = false;
                    }
                    else
                    {
                        // log all relevant data for this completed section
                        lohSVsMatchedCount += writeLOHData(sampleId, currentChr, lohStartCnData, cnData, priorCN, lohMinCN,
                                lohSegments, false, false, true);
                        ++lohSectionCount;
                    }

                    reset = true;
                }
                else if(cnData.matchesSegment(TELOMERE, false))
                {
                    // rest of arm was lost so no linking SV for LOH section - but still record the event
                    writeLOHData(sampleId, currentChr, lohStartCnData, cnData, priorCN, lohMinCN,
                            lohSegments, true, false, true);
                    reset = true;
                }
                else if(cnData.CopyNumber < TOTAL_CN_LOSS)
                {
                    // other chromatid loss has occurred - this will cancel a valid LOH due to uncertainty unless due to a simple SV
                    if(cnData.matchesSegment(SegmentSupport.DEL, true) && cnData.matchesSegment(SegmentSupport.DEL, false)
                    && isSingleVariant(cnData))
                    {
                        LOGGER.debug("total CN loss matches single SV({} : {} -> {})",
                                currentChr, cnData.StartPos, cnData.EndPos);
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
            else if(!reset)
            {
                boolean lohLost = minCN < MIN_LOH_CN;

                // check that an SV with correct orientation exists here
                if(lohLost && cnData.matchesSV(true) && findSvData(cnData, 1) == null)
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

                    if(lohOnStartTelomere)
                    {
                        LOGGER.debug("chr({}) LOH at telemore", currentChr);
                    }
                    else
                    {
                        LOGGER.debug(String.format("chr(%s) starting LOH at pos(%d) minCN(%.3f cn=%.3f baf=%.3f) priorCN(%.3f)",
                                currentChr, cnData.StartPos, cnData.CopyNumber, cnData.ActualBaf, lohMinCN, priorCN));

                        // check for segments ending on telomere
                        if(cnData.matchesSegment(TELOMERE, false))
                        {
                            writeLOHData(sampleId, currentChr, lohStartCnData, lohStartCnData, priorCN, lohMinCN,
                                    lohSegments, true, false, true);
                            reset = true;
                        }
                    }
                }
            }

            if(reset)
            {
                isLohSection = false;
                lohOnStartTelomere = false;
                lohSegments = 0;
                lohStartCnData = null;
                totalLoss = false;

                if(newChromosome && (minCN < MIN_LOH_CN))
                {
                    lohOnStartTelomere = true;
                    lohStartCnData = cnData;
                }

                currentChr = cnData.Chromosome;
            }

            lastMinCN = minCN;
        }

        LOGGER.info("sample({}) LOH sections({}) fullMatched({})",
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

    private StructuralVariantData findSvData_old(final SvCNData cnData, int requiredOrient)
    {
        if (cnData.matchesSegment(UNKNOWN, true)
                || cnData.matchesSegment(TELOMERE, true)
                || cnData.matchesSegment(CENTROMERE, true))
        {
            return null;
        }

        long svPosition = requiredOrient == -1 ? cnData.StartPos : cnData.StartPos - 1;

        for(final StructuralVariantData var : mSvDataList)
        {
            if(var.filter().equals(PON_FILTER_PON))
                continue;

            if(var.startChromosome().equals(cnData.Chromosome) && var.startOrientation() == requiredOrient && var.startPosition() == svPosition)
            {
                return var;
            }

            if(var.endChromosome().equals(cnData.Chromosome) && var.endOrientation() == requiredOrient && var.endPosition() == svPosition)
            {
                return var;
            }
        }

        return null;
    }

    // private static int SPEC_CN_ID = -1;
    private static int SPEC_CN_ID = 2317876;

    private int writeLOHData(final String sampleId, final String chr, SvCNData startData, SvCNData endData,
            double lastMinCN, double lohMinCN, int segCount, boolean incomplete, boolean skipped, boolean isValid)
    {
        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputPath;

                outputFileName += "CN_LOH_EVENTS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                // SV info
                mFileWriter.write("SampleId,Chromosome,CnIdStart,CnIdEnd,PosStart,PosEnd,SegStart,SegEnd,");
                mFileWriter.write("PrevCN,StartCN,EndCN,MinCN,SegCount,Length,StartSV,EndSV,Skipped,IsValid");
                mFileWriter.newLine();
            }

            StructuralVariantData startSvData = findSvData(startData, !skipped ? 1 : -1);

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
                endSvData = findSvData(endData, !skipped ? -1 : 1);
            }
            else
            {
                // segment has either started and/or finished on the telomere segment
                endData = startData;
                lohLength = startData.EndPos - startData.StartPos + 1;
            }

            if(lohLength <= 0)
            {
                LOGGER.warn("negative length({})", lohLength);
            }

            if (startSvData != null && endSvData != null)
            {
                if (startSvData.id().equals(endSvData.id()))
                {
                    LOGGER.debug("sample({}) cnID({} -> {}) matches singleSV({} - {})",
                            sampleId, startData.asString(), endData.asString(), startSvData.id(), startSvData.type());
                }
                else
                {
                    LOGGER.debug("sample({}) cnID({} -> {}) matches pairSV({} -> {})",
                            sampleId, startData.asString(), endData.asString(), startSvData.id(), endSvData.id());
                }
            }
            else
            {
                if(!skipped && !incomplete && isValid)
                {
                    if (startSvData == null && startData.matchesSV(true))
                    {
                        boolean incorrectOrientation = (startData.getStructuralVariantData() != null);

                        LOGGER.debug(String.format("sample(%s) LOH start segment(%s) no SV match: orient(%s) type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                                sampleId, startData.asString(), incorrectOrientation ? "wrong" : "ok", startData.SegStart,
                                startData.ActualBaf, startData.ObservedBaf, startData.BafCount, startData.CopyNumber));
                    }

                    if (!incomplete && endSvData == null && endData.matchesSV(true))
                    {
                        boolean incorrectOrientation = (endData.getStructuralVariantData() != null);

                        LOGGER.debug(String.format("sample(%s) LOH end segment(%s) no SV match: orient(%s), type(%s) baf(actual=%.4f observed=%.4f count=%d) copyNumber(%.4f)",
                                sampleId, endData.asString(), incorrectOrientation ? "wrong" : "ok", endData.SegStart,
                                endData.ActualBaf, endData.ObservedBaf, endData.BafCount, endData.CopyNumber));
                    }
                }

                LOGGER.debug("sample({}) cnID({} -  -> {}) not fully matched pairSV({} -> {})",
                        sampleId, startData.asString(), startData.SegStart, endData.asString(), endData.SegStart,
                        startSvData != null ? startSvData.id() : "", endSvData != null ? endSvData.id() : "");
            }

            mFileWriter.write(String.format("%s,%s,%d,%d,%d,%d,%s,%s",
                    sampleId, chr, startData.id(), endData.id(),
                    startData.StartPos, incomplete ? endData.EndPos : endData.StartPos,
                    startData.SegStart, incomplete ? endData.SegEnd : endData.SegStart));

            mFileWriter.write(String.format(",%.4f,%.4f,%.4f,%.4f,%d,%d,%s,%s,%s,%s",
                    lastMinCN, (1 - startData.ActualBaf) * startData.CopyNumber, (1 - endData.ActualBaf) * endData.CopyNumber,
                    lohMinCN, segCount, lohLength,
                    startSvData != null ? startSvData.id() : "0", endSvData != null ? endSvData.id() : "0",
                    skipped, isValid));

            mFileWriter.newLine();

            return (startSvData != null && endSvData != null) ? 1 : 0;
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to copy number LOH outputFile: {}", e.toString());
            return 0;
        }
    }

    private static int LOH_COLUMN_COUNT = 18;

    public void loadLOHFromCSV(final String filename, final String specificSample)
    {
        if (filename.isEmpty())
            return;

        mSampleLohData = new HashMap();

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty LOH CSV file({})", filename);
                return;
            }

            String currentSample = "";
            List<SvLOH> lohDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != LOH_COLUMN_COUNT)
                    continue;

                String sampleId = items[0];

                if(!specificSample.isEmpty() && !specificSample.equals(sampleId))
                {
                    continue;
                }

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    if(!currentSample.isEmpty())
                        mSampleLohData.put(currentSample, lohDataList);

                    lohDataList = Lists.newArrayList();
                    currentSample = sampleId;
                }

                // CSV fields
                // SampleId,Chromosome,CnIdStart,CnIdEnd,PosStart,PosEnd,SegStart,SegEnd,PrevCN,StartCN,EndCN,MinCN,SegCount,Length,StartSV,EndSV
                // 0  1        2   3        4      5        6      7        8           9         10         11

                SvLOH lohData = new SvLOH(
                        items[0],
                        items[1],
                        Integer.parseInt(items[2]),
                        Integer.parseInt(items[3]),
                        Long.parseLong(items[4]),
                        Long.parseLong(items[5]),
                        items[6],
                        items[7],
                        Double.parseDouble(items[8]),
                        Double.parseDouble(items[9]),
                        Double.parseDouble(items[10]),
                        Double.parseDouble(items[11]),
                        Integer.parseInt(items[12]),
                        Long.parseLong(items[13]),
                        items[14],
                        items[15],
                        Boolean.parseBoolean(items[16]),
                        Boolean.parseBoolean(items[17]));

                lohDataList.add(lohData);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read LOH CSV file({}): {}", filename, e.toString());
        }
    }

    public static int P_ARM_TELOMERE_CN = 0;
    public static int CENTROMERE_CN = 1;
    public static int Q_ARM_TELOMERE_CN = 2;

    public final Map<String, double[]> createChrCopyNumberMap(final List<PurpleCopyNumber> cnRecords)
    {
        final Map<String, double[]> chrMap = new HashMap();

        double[] currentChrData = null;

        for(final PurpleCopyNumber cnRecord : cnRecords)
        {
            if(cnRecord.segmentStartSupport().equals(TELOMERE))
            {
                currentChrData = new double[Q_ARM_TELOMERE_CN+1];
                chrMap.put(cnRecord.chromosome(), currentChrData);
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

        return chrMap;
    }

    public final List<StructuralVariantData> createNoneSegments(final List<PurpleCopyNumber> cnRecords)
    {
        List<StructuralVariantData> svList = Lists.newArrayList();

        for(int i = 0; i < cnRecords.size(); ++i)
        {
            if(i + 1 >= cnRecords.size())
                break;

            final PurpleCopyNumber prevCnRecord = cnRecords.get(i);
            final PurpleCopyNumber noneCnRecord = cnRecords.get(i+1);

            if(prevCnRecord.segmentEndSupport() != NONE || noneCnRecord.segmentStartSupport() != NONE)
                continue;

            double segmentCopyNumber = noneCnRecord.averageTumorCopyNumber();
            double copyNumberDiff = segmentCopyNumber - prevCnRecord.averageTumorCopyNumber();
            double copyNumberChange = abs(copyNumberDiff);

            double copyNumber = 0;

            if(copyNumberDiff > 0)
            {
                if(i < cnRecords.size() - 1)
                {
                    final PurpleCopyNumber nextCnRecord = cnRecords.get(i+1);
                    copyNumber = nextCnRecord.averageTumorCopyNumber();
                }
            }
            else
            {
                copyNumber = prevCnRecord.averageTumorCopyNumber();
            }

            // negative CN change means sequence has come in from left (a lower position) and dropped at the breakend
            byte orientation = copyNumberDiff > 0 ? (byte)-1 : 1;

            long position = orientation == -1 ? noneCnRecord.start() : noneCnRecord.start() - 1;

            int varId = mNoneSvId++;

            svList.add(
                    ImmutableStructuralVariantData.builder()
                            .id(Integer.toString(varId))
                            .vcfId("")
                            .type(SGL)
                            .ploidy(copyNumberChange)
                            .ploidyMin(copyNumberChange)
                            .ploidyMax(copyNumberChange)
                            .startPosition(position)
                            .startChromosome(noneCnRecord.chromosome())
                            .startOrientation(orientation)
                            .startAF(0.0)
                            .adjustedStartAF(0.0)
                            .adjustedStartCopyNumber(copyNumber)
                            .adjustedStartCopyNumberChange(copyNumberChange)
                            .homology("")
                            .inexactHomologyOffsetStart(0)
                            .insertSequence("")
                            .imprecise(false)
                            .qualityScore(0.0)
                            .filter(NONE_SEGMENT_INFERRED)
                            .endPosition(-1)
                            .endChromosome("0")
                            .endOrientation((byte)1)
                            .endAF(0.0)
                            .adjustedEndAF(0.0)
                            .adjustedEndCopyNumber(0.0)
                            .adjustedEndCopyNumberChange(0.0)
                            .event("")
                            .startTumourVariantFragmentCount(0)
                            .startTumourReferenceFragmentCount(0)
                            .startNormalVariantFragmentCount(0)
                            .startNormalReferenceFragmentCount(0)
                            .endTumourVariantFragmentCount(0)
                            .endTumourReferenceFragmentCount(0)
                            .endNormalVariantFragmentCount(0)
                            .endNormalReferenceFragmentCount(0)
                            .startIntervalOffsetStart(0)
                            .startIntervalOffsetEnd(0)
                            .endIntervalOffsetStart(0)
                            .endIntervalOffsetEnd(0)
                            .inexactHomologyOffsetStart(0)
                            .inexactHomologyOffsetEnd(0)
                            .startLinkedBy("")
                            .endLinkedBy("")
                            .startRefContext("")
                            .endRefContext("")
                            .insertSequenceAlignments("")
                            .recovered(false)
                            .build());
        }

        return svList;
    }

    private void initialisePloidyWriter()
    {
        try
        {
            if (mRecalcPloidyFileWriter == null)
            {
                String outputFileName = mOutputPath;

                outputFileName += "CN_PLOIDY_CALC_DATA.csv";

                mRecalcPloidyFileWriter = createBufferedWriter(outputFileName, false);

                if(!mWriteVerbosePloidyData)
                {
                    // SV info
                    mRecalcPloidyFileWriter.write("SampleId,SvId,EstPloidy,EstUncertainty");
                }
                else
                {
                    mRecalcPloidyFileWriter.write("SampleId,SvId,Type,Ploidy,TumorReadCount");
                    mRecalcPloidyFileWriter.write(",ChrStart,PosStart,OrientStart,CNStart,CNChgStart,DWCountStart,PrevDWCountStart,NextDWCountStart");
                    mRecalcPloidyFileWriter.write(",ChrEnd,PosEnd,OrientEnd,CNEnd,CNChgEnd,DWCountEnd,PrevDWCountEnd,NextDWCountEnd");
                    mRecalcPloidyFileWriter.write(",PoisRcLow,PoisRcHigh,PloidyUncertainty,EstPloidy,EstUncertainty,MinPloidy,MaxPloidy");
                }

                mRecalcPloidyFileWriter.newLine();
            }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to ploidy recalc outputFile: {}", e.toString());
        }
    }

    private void reaclcAdjustedPloidy(final String sampleId)
    {
        if (mWriteAdjustedPloidyToFile)
        {
            initialisePloidyWriter();
        }

        List<StructuralVariantData> updatedSvDataList = Lists.newArrayList();

        try
        {
            BufferedWriter writer = mRecalcPloidyFileWriter;

            for (Map.Entry<String, SvCNData[]> entry : mSvCnDataMap.entrySet())
            {
                final SvCNData[] cnDataPair = entry.getValue();

                final SvCNData cnStartData = cnDataPair[SVI_START];

                if (cnStartData == null)
                    continue;

                final SvCNData cnEndData = cnDataPair[SVI_END]; // may be null

                final StructuralVariantData svData = cnStartData.getStructuralVariantData();

                double ploidy = svData.ploidy();
                int tumorReadCount = svData.startTumourVariantFragmentCount();

                final int[] startDepthData = extractDepthWindowData(cnStartData);
                double cnStart = svData.adjustedStartCopyNumber();
                double cnChgStart = svData.adjustedStartCopyNumberChange();

                int[] endDepthData = null;
                double cnEnd = 0;
                double cnChgEnd = 0;

                if (cnEndData != null)
                {
                    endDepthData = extractDepthWindowData(cnEndData);
                    cnEnd = svData.adjustedEndCopyNumber();
                    cnChgEnd = svData.adjustedEndCopyNumberChange();
                }
                else
                {
                    endDepthData = new int[2];
                }

                final double calcResults[] = calcAdjustedPloidyValues(cnStart, cnChgStart, cnEnd, cnChgEnd,
                        ploidy, tumorReadCount, startDepthData, cnEndData != null ? endDepthData : null);

                double ploidyEstimate = calcResults[APC_EST_PLOIDY];
                double ploidyUncertainty = calcResults[APC_EST_UNCERTAINTY];

                // NONE segment SVs will usually fail on the ploidy recalc
                if(!svData.filter().equals(NONE_SEGMENT_INFERRED) && (Double.isNaN(ploidyEstimate) || Double.isNaN(ploidyUncertainty)))
                {
                    LOGGER.debug("sample({}) svID({} type={}) ploidy(est={} unc={}",
                            sampleId, svData.id(), svData.type(), ploidyEstimate, ploidyUncertainty);
                }

                if(Double.isNaN(ploidyEstimate))
                    ploidyEstimate = 0;

                if(Double.isNaN(ploidyUncertainty))
                    ploidyUncertainty = 0;

                updatedSvDataList.add(ImmutableStructuralVariantData.builder()
                        .from(svData)
                        .ploidyMin(ploidyEstimate - ploidyUncertainty)
                        .ploidyMax(ploidyEstimate + ploidyUncertainty)
                        .build());

                if (writer != null)
                {
                    if (!mWriteVerbosePloidyData)
                    {
                        writer.write(String.format("%s,%s,%.4f,%.4f",
                                sampleId, svData.id(),
                                ploidyEstimate, ploidyUncertainty));
                    }
                    else
                    {
                        writer.write(String.format("%s,%s,%s,%.4f,%d",
                                sampleId, svData.id(), svData.type(), svData.ploidy(), tumorReadCount));

                        writer.write(String.format(",%s,%d,%d,%.4f,%.4f,%d,%d,%d",
                                svData.startChromosome(), svData.startPosition(), svData.startOrientation(),
                                cnStart, cnChgStart, cnStartData.DepthWindowCount, startDepthData[0], startDepthData[1]));

                        writer.write(String.format(",%s,%d,%d,%.4f,%.4f,%d,%d,%d",
                                svData.endChromosome(), svData.endPosition(), svData.endOrientation(),
                                cnEnd, cnChgEnd, cnEndData != null ? cnEndData.DepthWindowCount : 0,
                                endDepthData != null ? endDepthData[0] : 0, endDepthData != null ? endDepthData[1] : 0));

                        writer.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                                calcResults[APC_RC_POIS_LOW], calcResults[APC_RC_POIS_HIGH],
                                calcResults[APC_PLOIDY_UNCERTAINTY],
                                ploidyEstimate, ploidyUncertainty,
                                ploidyEstimate - ploidyUncertainty,
                                ploidyEstimate + ploidyUncertainty));
                    }

                    writer.newLine();
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to ploidy recalc outputFile: {}", e.toString());
        }

        if (mUpdateAdjustedPloidyToDB)
        {
            mDbAccess.updateCalculatedPloidy(updatedSvDataList);
        }
    }

    public static int APC_RC_POIS_LOW = 0;
    public static int APC_RC_POIS_HIGH = 1;
    public static int APC_PLOIDY_UNCERTAINTY = 2;
    public static int APC_EST_PLOIDY = 3;
    public static int APC_EST_UNCERTAINTY = 4;
    public static int APC_RESULTS_MAX = 5;

    private static double POIS_PROB_LOW = 0.005;
    private static double POIS_PROB_HIGH = 0.995;

    public static double[] calcAdjustedPloidyValues(double cnStart, double cnChgStart, double cnEnd, double cnChgEnd,
            double ploidy, int tumorReadCount, final int[] startDepthData, final int[] endDepthData)
    {
        double cnUncertaintyStart = calcCopyNumberSideUncertainty(cnStart, startDepthData);
        double cnUncertaintyEnd = endDepthData != null ? calcCopyNumberSideUncertainty(cnEnd, endDepthData) : 0;

        double poissonRCLow = calcPoisonReadCount(tumorReadCount, POIS_PROB_LOW);
        double poissonRCHigh = calcPoisonReadCount(tumorReadCount, POIS_PROB_HIGH);

        double ploidyUncertainty = ploidy * (poissonRCHigh - poissonRCLow) * 0.5 / tumorReadCount;

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

        observations.add(ploidy);
        uncertainties.add(ploidyUncertainty);

        double estPloidy;
        double estUncertainty;

        if(observations.size() == 1)
        {
            estPloidy = ploidy;
            estUncertainty = ploidyUncertainty;
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

        results[APC_RC_POIS_LOW] = poissonRCLow;
        results[APC_RC_POIS_HIGH] = poissonRCHigh;
        results[APC_PLOIDY_UNCERTAINTY] = ploidyUncertainty;
        results[APC_EST_UNCERTAINTY] = estUncertainty;
        results[APC_EST_PLOIDY] = estPloidy;

        return results;
    }

    private static double ABS_UNCERTAINTY = 0.15;
    private static double RELATIVE_UNCERTAINTY = 0.1;
    private static double ADDITIONAL_ABS_UNCERTAINTY = 0.4;
    private static double ADDITIONAL_REL_UNCERTAINTY = 0.15;

    private static double calcCopyNumberSideUncertainty(double copyNumber, final int[] depthData)
    {
        int minDepthCount = min(depthData[0], depthData[1]);

        if(minDepthCount <= 0)
            return 0;

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

        /*
        double initRange = 3.7 / sqrt(expectedVal); // works for requiredProb = 1e-4
        int testValue = (int) Double.max(round(expectedVal * (1 - initRange)), 0);
        int testValueUpper = (int) Double.max(round(expectedVal * (1 - initRange*0.5)), 0);
        int testValueLower = (int) Double.max(round(expectedVal * (1 - initRange*2)), 0);
        */

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


    private int[] extractDepthWindowData(final SvCNData cnData)
    {
        // get depth window count before, at and after this segment
        int[] deptWindowData = new int[2];

        if(cnData.getIndex() > 0)
        {
            final SvCNData prevCnData = mCnDataList.get(cnData.getIndex()-1);

            if(prevCnData.Chromosome.equals(cnData.Chromosome))
                deptWindowData[0] = prevCnData.DepthWindowCount;
        }

        if(cnData.getIndex() < mCnDataList.size() - 1)
        {
            final SvCNData nextCnData = mCnDataList.get(cnData.getIndex() + 1);

            if(nextCnData.Chromosome.equals(cnData.Chromosome))
                deptWindowData[1] = nextCnData.DepthWindowCount;
        }

        return deptWindowData;
    }

    private static int PLOIDY_CALC_COLUMN_COUNT = 4;

    public void loadPloidyCalcData(final String filename, final String specificSample)
    {
        if (filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty LOH CSV file({})", filename);
                return;
            }

            String currentSample = "";
            Map<String,double[]> svDataMap = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != PLOIDY_CALC_COLUMN_COUNT)
                    continue;

                String sampleId = items[0];

                if(!specificSample.isEmpty() && !specificSample.equals(sampleId))
                {
                    if(mSampleSvPloidyCalcMap.size() == 1)
                        break;

                    continue;
                }

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    svDataMap = new HashMap();
                    mSampleSvPloidyCalcMap.put(currentSample, svDataMap);
                }

                final String svId = items[1];
                double[] ploidyCalcs = {Double.parseDouble(items[2]), Double.parseDouble(items[3])};

                svDataMap.put(svId, ploidyCalcs);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read SV ploidy calcs file({}): {}", filename, e.toString());
        }
    }


    public void close()
    {
        closeBufferedWriter(mFileWriter);
        closeBufferedWriter(mRecalcPloidyFileWriter);
    }

}
