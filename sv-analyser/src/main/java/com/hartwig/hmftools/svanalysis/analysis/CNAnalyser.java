package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.svanalysis.SvAnalyser.SAMPLE;
import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PON;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_UNKNOWN;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_TELOMERE;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_CENTROMERE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
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
import org.jetbrains.annotations.NotNull;

public class CNAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(CNAnalyser.class);

    private final String mOutputPath;
    private final SvUtilities mUtils;
    private Map<String, List<SvCNData>> mSampleCNData;

    private List<String> mSampleIds;
    private int mRecordId;

    BufferedWriter mFileWriter;
    DatabaseAccess mDbAccess;
    List<StructuralVariantData> mSvDataList;

    Map<String, List<SvLOH>> mSampleLohData;

    private boolean mAnalyseLOH;
    private boolean mAnalyseFlips;

    private static final String COPY_NUMBER_FILE = "cn_file";

    private static double CN_ROUNDING= 0.2;
    private static double CN_DIFF_MARGIN = 0.25;
    private static double CN_CHANGE_MIN = 0.8;
    private static int DB_MAX_LENGTH = 1000;

    private static double MIN_LOH_CN = 0.5;

    public CNAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mSampleCNData = new HashMap();
        mUtils = new SvUtilities(0);
        mOutputPath = outputPath;
        mSampleIds = Lists.newArrayList();
        mFileWriter = null;
        mDbAccess = dbAccess;
        mSvDataList = Lists.newArrayList();
        mSampleLohData = null;
        mRecordId = 0;

        mAnalyseLOH = true;
        mAnalyseFlips = false;
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

        setAnalyseFlips(false);
        setAnalyseLOH(true);

        return true;
    }

    public void setAnalyseLOH(boolean toggle) { mAnalyseLOH = toggle; }
    public void setAnalyseFlips(boolean toggle) { mAnalyseFlips = toggle; }

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

        if(mAnalyseFlips)
            analyseFlips(sampleId, sampleData, "");
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
        if(cnData.segStart().equals(CN_SEG_NONE) || cnData.segStart().equals(CN_SEG_UNKNOWN)
        || cnData.segStart().equals(CN_SEG_TELOMERE) || cnData.segStart().equals(CN_SEG_CENTROMERE))
            return null;

        for(final StructuralVariantData var : mSvDataList)
        {
            if(var.filter().equals(PON_FILTER_PON))
                continue;

            // if(svData.type() == StructuralVariantType.BND)
            if(var.startChromosome().equals(cnData.chromosome()) && var.startPosition() == cnData.startPos())
            {
                if(requiredOrient == 0 || var.startOrientation() == requiredOrient)
                    return var;
            }

            if(var.endChromosome().equals(cnData.chromosome()) && var.endPosition() == cnData.startPos())
            {
                if(requiredOrient == 0 || var.endOrientation() == requiredOrient)
                    return var;
            }
        }

        return null;
    }

    public final Map<String, List<SvLOH>> getSampleLohData() { return mSampleLohData; }

    private static int LOH_COLUMN_COUNT = 16;

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

            int cnCount = 0;
            String currentSample = "";
            List<SvLOH> lohDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != LOH_COLUMN_COUNT) {
                    continue;
                }

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
                        items[15]);

                lohDataList.add(lohData);
            }
        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read LOH CSV file({})", filename);
        }
    }

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
        double priorCN = 0; // before the LOGH segment
        int lohSegments = 0;
        SvCNData lohStartCN = null;
        int lohSectionCount = 0;
        int lohSVsMatchedCount = 0;

        for(SvCNData cnData : cnDataList)
        {
            double minCN = (1 - cnData.actualBaf()) * cnData.copyNumber();

            boolean reset = currentChr.isEmpty() || (!currentChr.isEmpty() && !cnData.chromosome().equals(currentChr));

            if(isLohSection)
            {
                if(minCN >= MIN_LOH_CN || cnData.segEnd().equals(CN_SEG_TELOMERE) || reset)
                {
                    // log all relevant data for this completed section
                    lohSVsMatchedCount += writeLOHData(sampleId, currentChr, lohStartCN, cnData, priorCN, lohMinCN, lohSegments);
                    ++lohSectionCount;

                    reset = true;
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

                    LOGGER.debug(String.format("chr(%s) starting LOH at pos(%d) minCN(%.3f cn=%.3f baf=%.3f) priorCN(%.3f)",
                            currentChr, cnData.startPos(), cnData.copyNumber(), cnData.actualBaf(), lohMinCN, priorCN));
                }
            }

            if(reset)
            {
                isLohSection = false;
                lohSegments = 0;
                lohStartCN = null;

                currentChr = cnData.chromosome();
                //currentArm = SvUtilities.CHROMOSOME_ARM_P;
            }

            lastMinCN = minCN;
        }

        LOGGER.info("sample({}) LOH sections({}) fullMatched({})",
                sampleId, lohSectionCount, lohSVsMatchedCount);
    }

    private int writeLOHData(
            final String sampleId, final String chr, SvCNData startData, SvCNData endData, double lastMinCN, double lohMinCN, int segCount)
    {
        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputPath;

                if (!outputFileName.endsWith("/"))
                    outputFileName += "/";

                outputFileName += "CN_LOH_ANALYSIS.csv";

                Path outputFile = Paths.get(outputFileName);

                mFileWriter = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

                // SV info
                mFileWriter.write("SampleId,Chromosome,CnIdStart,CnIdEnd,PosStart,PosEnd,SegStart,SegEnd,");
                mFileWriter.write("PrevCN,StartCN,EndCN,MinCN,SegCount,Length,StartSV,EndSV");
                mFileWriter.newLine();
            }

            StructuralVariantData startSvData = findSvData(startData, 1);

            StructuralVariantData endSvData = null;
            long lohLength = 0;

            if(endData.chromosome().equals(chr) && !startData.segEnd().equals(CN_SEG_TELOMERE))
            {
                lohLength = endData.startPos() - startData.startPos();
                endSvData = findSvData(endData, -1);
            }
            else
            {
                // segment has either started and finished on the last (telomere) segement or finished the next chromosome
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
                    startData.segStart(), endData.segStart()));

            mFileWriter.write(String.format(",%.4f,%.4f,%.4f,%.4f,%d,%d,%s,%s",
                    lastMinCN, (1 - startData.actualBaf()) * startData.copyNumber(), (1 - endData.actualBaf()) * endData.copyNumber(),
                    lohMinCN, segCount, lohLength,
                    startSvData != null ? startSvData.id() : "0", endSvData != null ? endSvData.id() : "0"));

            mFileWriter.newLine();

            return (startSvData != null && endSvData != null) ? 1 : 0;
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
            return 0;
        }
    }

    private void analyseFlips(final String sampleId, List<SvCNData> cnDataList, final String specificChromosome)
    {
        int segCount = 0;
        double lastCNChange = 0;
        double startCN = 0;
        double maxCN = 0;
        double minCN = 0;
        int delCount = 0;
        int dupCount = 0;
        int invCount = 0;
        int bndCount = 0;
        String currentChr = "";
        String currentArm = SvUtilities.CHROMOSOME_ARM_P;
        boolean checkCNChange = false;
        List<Long> dbLengths = Lists.newArrayList();
        Map<Double, Integer> cnChangeMap = new HashMap();

        if(mDbAccess != null)
            loadSVData(sampleId);

        SvCNData lastCnData = null;

        for(SvCNData cnData : cnDataList)
        {
            if(!specificChromosome.isEmpty() && !specificChromosome.equals(cnData.chromosome()))
            {
                if(segCount > 0)
                    break;

                continue;
            }

            double copyNumber = cnData.copyNumber();

            if(currentChr.isEmpty()
            || (!currentChr.isEmpty() && !cnData.chromosome().equals(currentChr))
            || (!currentArm.isEmpty() && cnData.segStart().equals(CN_SEG_CENTROMERE)))
            {
                if(!currentChr.isEmpty())
                {
                    writeFlipData(
                            sampleId, currentChr, currentArm, startCN, minCN, maxCN, segCount, dbLengths, cnChangeMap,
                            delCount, dupCount, invCount, bndCount);
                }

                // reset as a new chromosome has been found
                segCount = 0;
                lastCNChange = 0;
                maxCN = copyNumber;
                minCN = copyNumber;
                startCN = copyNumber;
                delCount = 0;
                dupCount = 0;
                invCount = 0;
                bndCount = 0;
                dbLengths.clear();
                cnChangeMap.clear();

                if(!currentChr.equals(cnData.chromosome()))
                {
                    currentChr = cnData.chromosome();
                    currentArm = SvUtilities.CHROMOSOME_ARM_P;
                }
                else
                {
                    currentArm = SvUtilities.CHROMOSOME_ARM_Q;
                }

                checkCNChange = false;
            }
            else
            {
                minCN = min(copyNumber, minCN);
                maxCN = max(copyNumber, maxCN);

                // analyse the change
                double thisCNChange = copyNumber - lastCnData.copyNumber();

                if(checkCNChange)
                {

                    if (lastCNChange < 0 && thisCNChange > 0
                        && abs(thisCNChange) >= CN_CHANGE_MIN && abs(thisCNChange + lastCNChange) <= CN_DIFF_MARGIN)
                    {
                        LOGGER.debug("sample({}) cnID({} -> {}) flipped cnChange({} -> {})", sampleId, lastCnData.asString(), cnData.asString(), lastCNChange, thisCNChange);
                        double cnRounded = roundCopyNumber(lastCNChange);

                        if (!cnChangeMap.containsKey(cnRounded))
                        {
                            cnChangeMap.put(cnRounded, 1);
                        }
                        else
                        {
                            cnChangeMap.replace(cnRounded, cnChangeMap.get(cnRounded) + 1);
                        }

                        long gapLength = cnData.startPos() - lastCnData.startPos();

                        // a deletion bridge (DB) is a CN drop at the last segment and regain on this segment, resulting from 2 distinct SVs
                        // so take the length of the last segment
                        if(!lastCnData.segStart().equals(StructuralVariantType.BND.toString()))
                        {
                            StructuralVariantData startSvData = findSvData(lastCnData, 1);
                            StructuralVariantData endSvData = findSvData(cnData, -1);

                            if (startSvData != null && endSvData != null) {

                                if (startSvData.id().equals(endSvData.id())) {
                                    LOGGER.debug("sample({}) cnID({} -> {}) matches singleSV({} - {})",
                                            sampleId, lastCnData.asString(), cnData.asString(), startSvData.id(), startSvData.type());
                                } else {

                                    LOGGER.debug("sample({}) cnID({} -> {}) matches pairSV({} -> {})",
                                            sampleId, lastCnData.asString(), cnData.asString(), startSvData.id(), endSvData.id());

                                    LOGGER.debug("sample({}) cnID({} -> {}) DB length({})",
                                            sampleId, lastCnData.asString(), cnData.asString(), gapLength);

                                    dbLengths.add(gapLength);

                                }
                            }
                            else
                            {
                                LOGGER.debug("sample({}) cnID({} -  -> {}) not fully matched pairSV({} -> {})",
                                        sampleId, lastCnData.asString(), lastCnData.segStart(), cnData.asString(), cnData.segStart(),
                                        startSvData != null ? startSvData.id() : "", endSvData != null ? endSvData.id() : "");
                            }
                        }

                        checkCNChange = false;
                    }
                }
                else
                {
                    checkCNChange = true;
                }

                lastCNChange = thisCNChange;
            }

            if(cnData.segStart().equals(StructuralVariantType.DEL.toString()))
                ++delCount;
            else if(cnData.segStart().equals(StructuralVariantType.DUP.toString()))
                ++dupCount;
            else if(cnData.segStart().equals(StructuralVariantType.INV.toString()))
                ++invCount;
            else if(cnData.segStart().equals(StructuralVariantType.BND.toString()))
                ++bndCount;

            lastCnData = cnData;
            ++segCount;

        }

        // write final chromosomal data
        if(!currentChr.isEmpty())
        {
            writeFlipData(
                    sampleId, currentChr, currentArm, startCN, minCN, maxCN, segCount, dbLengths, cnChangeMap,
                    delCount, dupCount,  invCount, bndCount);
        }
    }
    
    private void writeFlipData(
            final String sampleId, final String chr, final String arm, double startCN, double minCN, double maxCN, int segCount,
            final List<Long> dbLengths, final Map<Double, Integer> cnChangeMap,
            int delCount, int dupCount, int invCount, int bndCount)
    {
        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputPath;

                if (!outputFileName.endsWith("/"))
                    outputFileName += "/";

                outputFileName += "CN_ANALYSIS.csv";

                Path outputFile = Paths.get(outputFileName);

                mFileWriter = Files.newBufferedWriter(outputFile); // , StandardOpenOption.CREATE_NEW

                // SV info
                mFileWriter.write("SampleId,Chromosome,Arm,StartCN,MinCN,MaxCN,SegCount,DelCount,DupCount,InvCount,BndCount,ArmLenRatio,");
                mFileWriter.write("FlipCount,MaxFlips,FlipValues,DBCount,ShortCount,DelLens");
                mFileWriter.newLine();
            }

            double armLengthRatio = 30000000.0 / mUtils.getChromosomalArmLength(chr, arm);

            mFileWriter.write(String.format("%s,%s,%s,%.2f,%.2f,%.2f,%d,%d,%d,%d,%d,%.2f",
                    sampleId, chr, arm, startCN, minCN, maxCN, segCount, delCount, dupCount, invCount, bndCount, armLengthRatio));

            int maxFlips = 0;
            double cnChgMax = 0;
            int totalFlips = 0;
            String cnChangesStr = "";

            for(Map.Entry<Double, Integer> entry : cnChangeMap.entrySet())
            {
                if(entry.getValue() > maxFlips)
                {
                    cnChgMax = entry.getKey();
                    maxFlips = entry.getValue();
                }

                totalFlips += entry.getValue();

                if(!cnChangesStr.isEmpty())
                    cnChangesStr += ";";

                cnChangesStr += String.format("%.1f=%d", entry.getKey(), entry.getValue());
            }

            // now include any other CN change values within the margin of error of the max value
            for(Map.Entry<Double, Integer> entry : cnChangeMap.entrySet())
            {
                double cnChange = entry.getKey();
                if (cnChange != cnChgMax && abs(cnChange - cnChgMax) <= CN_DIFF_MARGIN)
                {
                    maxFlips += entry.getValue();
                }
            }

            int shortDBCount = 0;
            String delLengthsStr = "";

            for(long delLen : dbLengths)
            {
                if(delLen <= DB_MAX_LENGTH)
                    ++shortDBCount;

                    if(!delLengthsStr.isEmpty())
                    delLengthsStr += ";";

                delLengthsStr += String.valueOf(delLen);

            }

            mFileWriter.write(String.format(",%d,%d,%s,%d,%d,%s",
                    totalFlips, maxFlips, cnChangesStr, dbLengths.size(), shortDBCount, delLengthsStr));

            mFileWriter.newLine();
        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
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
