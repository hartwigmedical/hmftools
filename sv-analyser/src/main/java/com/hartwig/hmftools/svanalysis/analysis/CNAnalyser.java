package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PON;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArmLength;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_UNKNOWN;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_TELOMERE;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_CENTROMERE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.amber.qc.AmberQCStatus;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import java.nio.file.StandardOpenOption;
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
    private Map<String, List<SvCNData>> mSampleCNData;

    private List<String> mSampleIds;
    private int mRecordId;
    private int mNoneSvId;

    BufferedWriter mFileWriter;
    DatabaseAccess mDbAccess;
    List<StructuralVariantData> mSvDataList;

    Map<String, List<SvLOH>> mSampleLohData;

    private boolean mAnalyseLOH;

    private static final String COPY_NUMBER_FILE = "cn_file";

    private static double CN_ROUNDING= 0.2;
    public static double MIN_LOH_CN = 0.5;

    public CNAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mSampleCNData = new HashMap();
        mOutputPath = outputPath;
        mSampleIds = Lists.newArrayList();
        mFileWriter = null;
        mDbAccess = dbAccess;
        mSvDataList = Lists.newArrayList();
        mSampleLohData = null;
        mRecordId = 0;
        mNoneSvId = 0;

        mAnalyseLOH = true;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(COPY_NUMBER_FILE, true, "Copy number CSV file");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> sampleIds)
    {
        mSampleIds.addAll(sampleIds);

        if(cmd.hasOption(COPY_NUMBER_FILE))
        {
            loadFromCSV(cmd.getOptionValue(COPY_NUMBER_FILE), mSampleIds.size() == 1 ? mSampleIds.get(0) : "");
        }

        setAnalyseLOH(true);

        return true;
    }

    public void setAnalyseLOH(boolean toggle) { mAnalyseLOH = toggle; }

    private void loadFromCSV(final String filename, final String specificSample)
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
                LOGGER.error("Empty copy number CSV file({})", filename);
                return;
            }

            int cnCount = 0;
            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 12) {
                    continue;
                }

                String sampleId = items[1];

                if(!specificSample.isEmpty() && !specificSample.equals(sampleId))
                {
                    continue;
                }

                if (!mSampleCNData.containsKey(sampleId))
                {
                    List<SvCNData> newList = Lists.newArrayList();
                    mSampleCNData.put(sampleId, newList);
                }

                List<SvCNData> sampleData = mSampleCNData.get(sampleId);

                // CSV fields
                // Id,SampleId,Chr,PosStart,PosEnd,SegStart,SegEnd,BafCount,ObservedBaf,ActualBaf,CopyNumber,CNMethod
                // 0  1        2   3        4      5        6      7        8           9         10         11

                SvCNData cnData = new SvCNData(
                        Integer.parseInt(items[0]),
                        items[2],
                        Long.parseLong(items[3]),
                        Long.parseLong(items[4]),
                        items[5],
                        items[6],
                        Integer.parseInt(items[7]),
                        Double.parseDouble(items[8]),
                        Double.parseDouble(items[9]),
                        Double.parseDouble(items[10]),
                        items[10]);

                sampleData.add(cnData);
                ++cnCount;
            }

            LOGGER.info("loaded  samples({}) copy number data count({})", mSampleCNData.size(), cnCount);

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read copy number CSV file({})", filename);
        }
    }

    public void analyseData()
    {
        int sampleCount = 0;

        if(!mSampleCNData.isEmpty())
        {
            for (Map.Entry<String, List<SvCNData>> entry : mSampleCNData.entrySet())
            {
                if(mSampleCNData.size() > 1)
                {
                    ++sampleCount;
                    LOGGER.info("analysing sample({}) with {} CN entries, totalProcessed({})",
                            entry.getKey(), entry.getValue().size(), sampleCount);
                }

                analyseData(entry.getKey(), entry.getValue());
            }
        }
        else if(mDbAccess != null)
        {
            for(final String sampleId : mSampleIds)
            {
                List<SvCNData> sampleData = loadFromDatabase(sampleId);

                ++sampleCount;
                LOGGER.info("analysing sample({}) with {} CN entries, totalProcessed({})",
                        sampleId, sampleData.size(), sampleCount);

                analyseData(sampleId, sampleData);
            }
        }
    }

    private List<SvCNData> loadFromDatabase(final String sampleId)
    {
        List<SvCNData> cnDataList = Lists.newArrayList();
        List<PurpleCopyNumber> cnRecords = mDbAccess.readCopynumbers(sampleId);

        for(final PurpleCopyNumber cnRecord : cnRecords)
        {
            cnDataList.add(new SvCNData(cnRecord, ++mRecordId));
        }

        return cnDataList;
    }

    private void analyseData(final String sampleId, final List<SvCNData> sampleData)
    {
        if(mAnalyseLOH)
            analyseLOH(sampleId, sampleData);
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

    public final List<StructuralVariantData> loadNoneSegments(final String sampleId)
    {
        List<StructuralVariantData> svList = Lists.newArrayList();

        List<PurpleCopyNumber> cnRecords = mDbAccess.readCopyNumberNoneSegments(sampleId);

        for(int i = 0; i < cnRecords.size(); ++i)
        {
            if(i + 1 >= cnRecords.size())
                break;

            final PurpleCopyNumber prevCnRecord = cnRecords.get(i);
            final PurpleCopyNumber noneCnRecord = cnRecords.get(i+1);

            if(prevCnRecord.segmentEndSupport() != SegmentSupport.NONE || noneCnRecord.segmentStartSupport() != SegmentSupport.NONE)
                continue;

            double copyNumber = noneCnRecord.averageTumorCopyNumber();
            double copyNumberDiff = copyNumber - prevCnRecord.averageTumorCopyNumber();
            double copyNumberChange = abs(copyNumberDiff);

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
                            .filter(AmberQCStatus.PASS.toString())
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

    private StructuralVariantData findSvData(final SvCNData cnData, int requiredOrient)
    {
        if(cnData.segStart().equals(CN_SEG_UNKNOWN)
        || cnData.segStart().equals(CN_SEG_TELOMERE) || cnData.segStart().equals(CN_SEG_CENTROMERE))
            return null;

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

    public final Map<String, List<SvLOH>> getSampleLohData() { return mSampleLohData; }

    private static int LOH_COLUMN_COUNT = 17;

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
                        Boolean.parseBoolean(items[16]));

                lohDataList.add(lohData);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read LOH CSV file({}): {}", filename, e.toString());
        }
    }

    // private static String SPECIFIC_CHR = "14";
    private static String SPECIFIC_CHR = "";
    private static int REMOTE_SV_DISTANCE = 1000000;

    private void analyseLOH(final String sampleId, List<SvCNData> cnDataList)
    {
        String currentChr = "";

        if(mDbAccess != null)
            loadSVData(sampleId);

        // walk through the CN records looking for any loss of hetrozygosity,
        // defined here as a change in the actual baf to zero
        // when CN rise back above zero, consider the section ended
        boolean isLohSection = false;
        double lohMinCN = 0;
        double lastMinCN = 0;
        double priorCN = 0; // before the LOH segment
        int lohSegments = 0;
        SvCNData lohStartCN = null;
        int lohSectionCount = 0;
        int lohSVsMatchedCount = 0;
        boolean lohOnStartTelomere = false;
        boolean totalLoss = false;

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

                            writeLOHData(sampleId, currentChr, cnData, nextData, priorCN, lohMinCN, 1, false, true);
                            continue;
                        }
                    }

                    if(lohOnStartTelomere || totalLoss)
                    {
                        // LOH section invalidated
                        if(lohOnStartTelomere)
                            writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN, lohSegments, false, false);

                        lohOnStartTelomere = false;
                        totalLoss = false;
                    }
                    else
                    {
                        // log all relevant data for this completed section
                        lohSVsMatchedCount += writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN, lohSegments, false, false);
                        ++lohSectionCount;
                    }

                    reset = true;
                }
                else if(cnData.segEnd().equals(CN_SEG_TELOMERE))
                {
                    // rest of arm was lost so no linking SV for LOH section - but still record the event
                    writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN, lohSegments, true, false);
                    reset = true;
                }
                else if(cnData.copyNumber() < 0.5)
                {
                    // other chromatid loss has occurred
                    totalLoss = true;
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
            double lastMinCN, double lohMinCN, int segCount, boolean incomplete, boolean skipped)
    {
        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputPath;

                if (!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "CN_LOH_EVENTS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                // SV info
                mFileWriter.write("SampleId,Chromosome,CnIdStart,CnIdEnd,PosStart,PosEnd,SegStart,SegEnd,");
                mFileWriter.write("PrevCN,StartCN,EndCN,MinCN,SegCount,Length,StartSV,EndSV,Skipped");
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
            else if(endData.chromosome().equals(chr) && !startData.segEnd().equals(CN_SEG_TELOMERE))
            {
                lohLength = endData.startPos() - startData.startPos();
                endSvData = findSvData(endData, !skipped ? -1 : 1);
            }
            else
            {
                // segment has either started and finished on the last (telomere) segment or finished the next chromosome
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
                    sampleId, chr, startData.id(), endData.id(), startData.startPos(), endData.startPos(),
                    startData.segStart(), incomplete ? endData.segEnd() : endData.segStart()));

            mFileWriter.write(String.format(",%.4f,%.4f,%.4f,%.4f,%d,%d,%s,%s,%s",
                    lastMinCN, (1 - startData.actualBaf()) * startData.copyNumber(), (1 - endData.actualBaf()) * endData.copyNumber(),
                    lohMinCN, segCount, lohLength,
                    startSvData != null ? startSvData.id() : "0", endSvData != null ? endSvData.id() : "0",
                    skipped));

            mFileWriter.newLine();

            return (startSvData != null && endSvData != null) ? 1 : 0;
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to copy number LOH outputFile: {}", e.toString());
            return 0;
        }
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

    private double roundCopyNumber(double copyNumber)
    {
        return abs(round(copyNumber / CN_ROUNDING) * CN_ROUNDING);
    }

}
