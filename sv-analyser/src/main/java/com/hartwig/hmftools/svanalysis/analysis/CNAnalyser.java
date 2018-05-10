package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.SvCNData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CNAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(CNAnalyser.class);

    private final String mOutputPath;
    private final SvUtilities mUtils;
    // private List<SvCNData> mAllCNData;
    private Map<String, List<SvCNData>> mSampleCNData;
    BufferedWriter mFileWriter;
    DatabaseAccess mDbAccess;
    List<StructuralVariantData> mSvDataList;

    private static double CN_ROUNDING= 0.2;
    private static double CN_DIFF_MARGIN = 0.25;
    private static double CN_CHANGE_MIN = 0.8;
    private static int DB_MAX_LENGTH = 1000;

    private static String TELOMORE = "TELOMERE";
    private static String CENTROMERE = "CENTROMERE";

    public CNAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mSampleCNData = new HashMap();
        mUtils = new SvUtilities(0);
        mOutputPath = outputPath;
        mFileWriter = null;
        mDbAccess = dbAccess;
        mSvDataList = Lists.newArrayList();
    }

    public void loadFromCSV(final String filename, final String specificSample)
    {
        if (filename.isEmpty()) {
            return;
        }

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null) {
                LOGGER.error("Empty copy number CSV file({})", filename);
                return;
            }

            // skip field names

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                // CSV fields
                // Id,SampleId,Chr,PosStart,PosEnd,SegStart,SegEnd,BafCount,Baf,CopyNumber,CNMethod
                // 0  1        2   3        4      5        6      7        8   9          10

                if (items.length != 11) {
                    continue;
                }

                String sampleId = items[1];

                if(!specificSample.isEmpty() && !specificSample.equals(sampleId))
                {
                    continue;
                }

                if (!mSampleCNData.containsKey(sampleId)) {
                    List<SvCNData> newList = Lists.newArrayList();
                    mSampleCNData.put(sampleId, newList);
                }

                List<SvCNData> sampleData = mSampleCNData.get(sampleId);

                SvCNData cnData = new SvCNData(
                        Integer.parseInt(items[0]),
                        items[2],
                        Long.parseLong(items[3]),
                        Long.parseLong(items[4]),
                        Double.parseDouble(items[9]),
                        items[5],
                        items[6],
                        items[10],
                        Double.parseDouble(items[8]));

                sampleData.add(cnData);
            }

            LOGGER.debug("loaded copy number data");

        } catch (IOException exception) {
            LOGGER.error("Failed to read copy number CSV file({})", filename);
        }
    }

    public void analyseData(final String specificSample, final String specificChromosome) {

        if(!specificSample.isEmpty())
        {
            if (!mSampleCNData.containsKey(specificSample)) {
                LOGGER.warn("sample({}) not found", specificSample);
                return;
            }

            List<SvCNData> sampleData = mSampleCNData.get(specificSample);

            analyseSample(specificSample, sampleData, specificChromosome);
        }
        else {
            for (Map.Entry<String, List<SvCNData>> entry : mSampleCNData.entrySet()) {

                LOGGER.info("analysing sample({}) with {} CN entries", entry.getKey(), entry.getValue().size());
                analyseSample(entry.getKey(), entry.getValue(), specificChromosome);

            }
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

    private StructuralVariantData findSvData(final SvCNData cnData)
    {
        if(cnData.getSegStart().equals("NONE"))
            return null;

        for(final StructuralVariantData svData : mSvDataList)
        {
            // if(svData.type() == StructuralVariantType.BND)
            if(svData.startChromosome().equals(cnData.getChromosome()) && svData.startPosition() == cnData.getStartPos())
            {
                return svData;
            }
            else if(svData.endChromosome().equals(cnData.getChromosome()) && svData.endPosition() == cnData.getStartPos())
            {
                return svData;
            }
        }

        return null;
    }

    private void analyseSample(final String sampleId, List<SvCNData> cnDataList, final String specificChromosome)
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
            if(!specificChromosome.isEmpty() && !specificChromosome.equals(cnData.getChromosome()))
            {
                if(segCount > 0)
                    break;

                continue;
            }

            double copyNumber = cnData.getCopyNumber();

            if(currentChr.isEmpty()
            || (!currentChr.isEmpty() && !cnData.getChromosome().equals(currentChr))
            || (!currentArm.isEmpty() && cnData.getSegStart().equals(CENTROMERE)))
            {
                if(!currentChr.isEmpty())
                {
                    writeSampleChromosomeData(
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

                if(!currentChr.equals(cnData.getChromosome()))
                {
                    currentChr = cnData.getChromosome();
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
                double thisCNChange = copyNumber - lastCnData.getCopyNumber();

                if(checkCNChange) {

                    if (lastCNChange < 0 && thisCNChange > 0
                        && abs(thisCNChange) >= CN_CHANGE_MIN && abs(thisCNChange + lastCNChange) <= CN_DIFF_MARGIN) {

                        LOGGER.debug("sample({}) cnID({} -> {}) flipped cnChange({} -> {})", sampleId, lastCnData.asString(), cnData.asString(), lastCNChange, thisCNChange);
                        double cnRounded = roundCopyNumber(lastCNChange);

                        if (!cnChangeMap.containsKey(cnRounded)) {
                            cnChangeMap.put(cnRounded, 1);
                        } else {
                            cnChangeMap.replace(cnRounded, cnChangeMap.get(cnRounded) + 1);
                        }

                        long gapLength = cnData.getStartPos() - lastCnData.getStartPos();

                        // a deletion bridge (DB) is a CN drop at the last segment and regain on this segment, resulting from 2 distinct SVs
                        // so take the length of the last segment
                        if(!lastCnData.getSegStart().equals(StructuralVariantType.BND.toString())) {

                            StructuralVariantData startSvData = findSvData(lastCnData);
                            StructuralVariantData endSvData = findSvData(cnData);

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
                            } else {
                                LOGGER.debug("sample({}) cnID({} -  -> {}) not fully matched pairSV({} -> {})",
                                        sampleId, lastCnData.asString(), lastCnData.getSegStart(), cnData.asString(), cnData.getSegStart(),
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

            if(cnData.getSegStart().equals(StructuralVariantType.DEL.toString()))
                ++delCount;
            else if(cnData.getSegStart().equals(StructuralVariantType.DUP.toString()))
                ++dupCount;
            else if(cnData.getSegStart().equals(StructuralVariantType.INV.toString()))
                ++invCount;
            else if(cnData.getSegStart().equals(StructuralVariantType.BND.toString()))
                ++bndCount;

            lastCnData = cnData;
            ++segCount;

        }

        // write final chromosomal data
        if(!currentChr.isEmpty())
        {
            writeSampleChromosomeData(
                    sampleId, currentChr, currentArm, startCN, minCN, maxCN, segCount, dbLengths, cnChangeMap,
                    delCount, dupCount,  invCount, bndCount);
        }

    }
    
    private void writeSampleChromosomeData(
            final String sampleId, final String chr, final String arm, double startCN, double minCN, double maxCN, int segCount,
            final List<Long> dbLengths, final Map<Double, Integer> cnChangeMap,
            int delCount, int dupCount, int invCount, int bndCount)
    {
        try {

            if (mFileWriter == null) {
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

            for(Map.Entry<Double, Integer> entry : cnChangeMap.entrySet()) {

                if(entry.getValue() > maxFlips) {
                    cnChgMax = entry.getKey();
                    maxFlips = entry.getValue();
                }

                totalFlips += entry.getValue();

                if(!cnChangesStr.isEmpty())
                    cnChangesStr += ";";

                cnChangesStr += String.format("%.1f=%d", entry.getKey(), entry.getValue());
            }

            // now include any other CN change values within the margin of error of the max value
            for(Map.Entry<Double, Integer> entry : cnChangeMap.entrySet()) {

                double cnChange = entry.getKey();
                if (cnChange != cnChgMax && abs(cnChange - cnChgMax) <= CN_DIFF_MARGIN) {
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

    public void close() {
        if (mFileWriter == null)
            return;

        try {
            mFileWriter.close();
        } catch (final IOException e) {
        }
    }

    private double roundCopyNumber(double copyNumber)
    {
        return abs(round(copyNumber / CN_ROUNDING) * CN_ROUNDING);
    }

}
