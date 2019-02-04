package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.NONE_SEGMENT_INFERRED;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvLOH;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CNAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(CNAnalyser.class);

    private final String mOutputPath;

    private List<String> mSampleIds;
    private int mRecordId;
    private int mNoneSvId;

    BufferedWriter mFileWriter;
    DatabaseAccess mDbAccess;
    List<StructuralVariantData> mSvDataList;

    Map<String, List<SvLOH>> mSampleLohData;

    private static final String COPY_NUMBER_FILE = "cn_file";

    public static double MIN_LOH_CN = 0.5;
    public static double TOTAL_CN_LOSS = 0.5;

    public CNAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mOutputPath = outputPath;
        mSampleIds = Lists.newArrayList();
        mFileWriter = null;
        mDbAccess = dbAccess;
        mSvDataList = Lists.newArrayList();
        mSampleLohData = null;
        mRecordId = 0;
        mNoneSvId = 0;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(COPY_NUMBER_FILE, true, "Copy number CSV file");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> sampleIds)
    {
        mSampleIds.addAll(sampleIds);
        return true;
    }

    public void findLOHEvents()
    {
        if(mDbAccess == null)
            return;

        int sampleCount = 0;

        for(final String sampleId : mSampleIds)
        {
            List<SvCNData> sampleData = loadCopyNumberData(sampleId);

            ++sampleCount;
            LOGGER.info("analysing sample({}) with {} CN entries, totalProcessed({})", sampleId, sampleData.size(), sampleCount);

            loadSVData(sampleId);
            findLohEvents(sampleId, sampleData);
        }
    }

    private List<SvCNData> loadCopyNumberData(final String sampleId)
    {
        List<SvCNData> cnDataList = Lists.newArrayList();
        List<PurpleCopyNumber> cnRecords = mDbAccess.readCopynumbers(sampleId);

        for(final PurpleCopyNumber cnRecord : cnRecords)
        {
            cnDataList.add(new SvCNData(cnRecord, ++mRecordId));
        }

        return cnDataList;
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

    private StructuralVariantData findSvData(final SvCNData cnData, int requiredOrient)
    {
        if(cnData.matchesSegment(UNKNOWN, true)
        || cnData.matchesSegment(TELOMERE, true)
        || cnData.matchesSegment(CENTROMERE, true))
        {
            return null;
        }

        long svPosition = requiredOrient == -1 ? cnData.startPos() : cnData.startPos() - 1;

        for(final StructuralVariantData var : mSvDataList)
        {
            if(var.filter().equals(PON_FILTER_PON))
                continue;

            if(var.startChromosome().equals(cnData.chromosome()) && var.startOrientation() == requiredOrient && var.startPosition() == svPosition)
            {
                return var;
            }

            if(var.endChromosome().equals(cnData.chromosome()) && var.endOrientation() == requiredOrient && var.endPosition() == svPosition)
            {
                return var;
            }
        }

        return null;
    }

    private boolean isSingleVariant(final SvCNData cnData)
    {
        for(final StructuralVariantData var : mSvDataList)
        {
            if(var.filter().equals(PON_FILTER_PON))
                continue;

            long svPosStart = var.startOrientation() == -1 ? cnData.startPos() : cnData.startPos() - 1;
            long svPosEnd = cnData.endPos() + 1;

            if(var.startChromosome().equals(cnData.chromosome()) && var.startPosition() == svPosStart && var.endPosition() == svPosEnd)
            {
                return true;
            }
        }

        return false;
    }

    public final Map<String, List<SvLOH>> getSampleLohData() { return mSampleLohData; }

    // private static String SPECIFIC_CHR = "13";
    private static String SPECIFIC_CHR = "";
    private static int REMOTE_SV_DISTANCE = 1000000;

    private void findLohEvents(final String sampleId, List<SvCNData> cnDataList)
    {
        String currentChr = "";

        // walk through the CN records looking for any loss of hetrozygosity, defined as a change in the actual baf to zero
        // and when CN rises back above this level, consider the section ended
        boolean isLohSection = false;
        double lohMinCN = 0;
        double lastMinCN = 0;
        double priorCN = 0; // before the LOH segment
        int lohSegments = 0;
        SvCNData lohStartCN = null;
        int lohSectionCount = 0;
        int lohSVsMatchedCount = 0;
        boolean lohOnStartTelomere = false;
        boolean totalLoss = false;StructuralVariantData totalLossSv = null;

        for(int index = 0; index < cnDataList.size(); ++index)
        {
            final SvCNData cnData = cnDataList.get(index);

            double minCN = (1 - cnData.actualBaf()) * cnData.copyNumber();

            boolean newChromosome = currentChr.isEmpty() || (!currentChr.isEmpty() && !cnData.chromosome().equals(currentChr));
            boolean reset = newChromosome;

            if(newChromosome && cnData.chromosome().equals(SPECIFIC_CHR))
            {
                LOGGER.debug("spec chr({})", SPECIFIC_CHR);
            }

            if(isLohSection || lohOnStartTelomere)
            {
                if(minCN >= MIN_LOH_CN || reset)
                {
                    // check for a short isolated TI and if found continue with the LOH
                    if(minCN >= MIN_LOH_CN && cnData.endPos() - cnData.startPos() <= SHORT_TI_LENGTH && index < cnDataList.size() - 1)
                    {
                        final SvCNData nextData = cnDataList.get(index+1);
                        double nextMinCN = (1 - nextData.actualBaf()) * nextData.copyNumber();

                        if(nextData.endPos() - cnData.startPos() > REMOTE_SV_DISTANCE && nextMinCN < MIN_LOH_CN
                        && lohStartCN!= null && cnData.startPos() - lohStartCN.startPos() > REMOTE_SV_DISTANCE)
                        {
                            LOGGER.debug("chr({}) skipping short isolated TI(id={} {} pos={} length={})",
                                    currentChr, cnData.id(), cnData.segStart(), cnData.startPos(), cnData.endPos() - cnData.startPos());

                            writeLOHData(sampleId, currentChr, cnData, nextData, priorCN, lohMinCN,
                                    1, false, true, false);
                            continue;
                        }
                    }

                    if(lohOnStartTelomere || totalLoss)
                    {
                        // LOH section invalidated
                        writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN,
                                lohSegments, false, false, !totalLoss);

                        lohOnStartTelomere = false;
                        totalLoss = false;
                    }
                    else
                    {
                        // log all relevant data for this completed section
                        lohSVsMatchedCount += writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN,
                                lohSegments, false, false, true);
                        ++lohSectionCount;
                    }

                    reset = true;
                }
                else if(cnData.matchesSegment(TELOMERE, false))
                {
                    // rest of arm was lost so no linking SV for LOH section - but still record the event
                    writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN,
                            lohSegments, true, false, true);
                    reset = true;
                }
                else if(cnData.copyNumber() < TOTAL_CN_LOSS)
                {
                    // other chromatid loss has occurred - this will cancel a valid LOH due to uncertainty unless due to a simple SV
                    if(cnData.matchesSegment(SegmentSupport.DEL, true)
                    && cnData.matchesSegment(SegmentSupport.DEL, false)
                    && isSingleVariant(cnData))
                    {
                        LOGGER.debug("total CN loss matches single SV({} : {} -> {})",
                                currentChr, cnData.startPos(), cnData.endPos());
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
                if (minCN < MIN_LOH_CN)
                {
                    // new LOH section identified
                    isLohSection = true;
                    lohSegments = 1;
                    lohStartCN = cnData;
                    lohMinCN = minCN;
                    priorCN = lastMinCN;

                    if(lohOnStartTelomere)
                    {
                        LOGGER.debug("chr({}) LOH at telemore", currentChr);
                    }
                    else
                    {
                        LOGGER.debug(String.format("chr(%s) starting LOH at pos(%d) minCN(%.3f cn=%.3f baf=%.3f) priorCN(%.3f)",
                                currentChr, cnData.startPos(), cnData.copyNumber(), cnData.actualBaf(), lohMinCN, priorCN));
                    }
                }
            }

            if(reset)
            {
                isLohSection = false;
                lohOnStartTelomere = false;
                lohSegments = 0;
                lohStartCN = null;
                totalLoss = false;

                if(newChromosome && (minCN < MIN_LOH_CN))
                {
                    lohOnStartTelomere = true;
                    lohStartCN = cnData;
                }

                currentChr = cnData.chromosome();
            }

            lastMinCN = minCN;
        }

        LOGGER.info("sample({}) LOH sections({}) fullMatched({})",
                sampleId, lohSectionCount, lohSVsMatchedCount);
    }

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
                lohLength = endData.endPos() - startData.startPos();
                endSvData = null;
            }
            else if(endData.chromosome().equals(chr) && !startData.matchesSegment(TELOMERE, false))
            {
                lohLength = endData.startPos() - startData.startPos();
                endSvData = findSvData(endData, !skipped ? -1 : 1);
            }
            else
            {
                // segment has either started and/or finished on the telomere segment
                endData = startData;
                lohLength = startData.endPos() - startData.startPos() + 1;
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
                LOGGER.debug("sample({}) cnID({} -  -> {}) not fully matched pairSV({} -> {})",
                        sampleId, startData.asString(), startData.segStart(), endData.asString(), endData.segStart(),
                        startSvData != null ? startSvData.id() : "", endSvData != null ? endSvData.id() : "");
            }

            mFileWriter.write(String.format("%s,%s,%d,%d,%d,%d,%s,%s",
                    sampleId, chr, startData.id(), endData.id(),
                    startData.startPos(), incomplete ? endData.endPos() : endData.startPos(),
                    startData.segStart(), incomplete ? endData.segEnd() : endData.segStart()));

            mFileWriter.write(String.format(",%.4f,%.4f,%.4f,%.4f,%d,%d,%s,%s,%s,%s",
                    lastMinCN, (1 - startData.actualBaf()) * startData.copyNumber(), (1 - endData.actualBaf()) * endData.copyNumber(),
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

            if(prevCnRecord.segmentEndSupport() != SegmentSupport.NONE || noneCnRecord.segmentStartSupport() != SegmentSupport.NONE)
                continue;

            double segmentCopyNumber = noneCnRecord.averageTumorCopyNumber();
            double copyNumberDiff = segmentCopyNumber - prevCnRecord.averageTumorCopyNumber();
            double copyNumberChange = abs(copyNumberDiff);
            double copyNumber = prevCnRecord.averageTumorCopyNumber();

            // negative CN change means sequence has come in from left (a lower position) and dropped at the breakend
            byte orientation = copyNumberDiff > 0 ? (byte)-1 : 1;

            long position = orientation == -1 ? noneCnRecord.start() : noneCnRecord.start() - 1;

            int varId = mNoneSvId++;

            svList.add(
                    ImmutableStructuralVariantData.builder()
                            .id(Integer.toString(varId))
                            .vcfId("")
                            .type(StructuralVariantType.SGL)
                            .ploidy(copyNumberChange)
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
                            .build());
        }

        return svList;
    }

    public void close()
    {
        if (mFileWriter == null)
            return;

        try
        {
            mFileWriter.close();
        }
        catch (final IOException e)
        {
        }
    }

}
