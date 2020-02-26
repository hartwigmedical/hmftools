package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.LONG;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

public class ExpectedRatesGenerator
{
    private final RnaExpConfig mConfig;

    // map of transcript (or unspliced) to all expected category counts covering it (and others)
    private final Map<String, List<TranscriptComboData>> mExpectedTransComboCounts;

    private ExpectedRatesData mCurrentExpRatesData;
    private int mCurrentFragSize;
    private int mCurrentFragFrequency;
    private int mCurrentReadLength;

    private final BufferedWriter mExpRateWriter;

    public static final String UNSPLICED_ID = "UNSPLICED";
    private static final String EXP_RATES_FILE_ID = "rna_trans_exp_rates";

    public static final int FL_SIZE = 0;
    public static final int FL_FREQUENCY = 1;

    public ExpectedRatesGenerator(final RnaExpConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mCurrentFragSize = 0;
        mCurrentFragFrequency = 0;
        mCurrentReadLength = mConfig.ReadLength;

        mExpectedTransComboCounts = Maps.newHashMap();
        mCurrentExpRatesData = null;

        mExpRateWriter = resultsWriter != null ? resultsWriter.getExpRatesWriter() : null;
    }

    public static ExpectedRatesGenerator from(final RnaExpConfig config)
    {
        return new ExpectedRatesGenerator(config, null);
    }

    public Map<String,List<TranscriptComboData>> getTransComboData() { return mExpectedTransComboCounts; }
    public ExpectedRatesData getExpectedRatesData() { return mCurrentExpRatesData; }

    public void generateExpectedRates(final GeneReadData geneReadData)
    {
        mExpectedTransComboCounts.clear();
        mCurrentExpRatesData = new ExpectedRatesData(geneReadData.GeneData.GeneId);

        final List<long[]> commonExonicRegions = geneReadData.getCommonExonicRegions();

        // apply fragment reads across each transcript as though it were fully transcribed
        final List<TranscriptData> transDataList = geneReadData.getTranscripts();

        for(final int[] flData : mConfig.ExpRateFragmentLengths)
        {
            mCurrentFragSize = flData[FL_SIZE];
            mCurrentFragFrequency = flData[FL_FREQUENCY];

            for (TranscriptData transData : geneReadData.getTranscripts())
            {
                boolean endOfTrans = false;

                for (ExonData exon : transData.exons())
                {
                    for (long startPos = exon.ExonStart; startPos <= exon.ExonEnd; ++startPos)
                    {
                        if (!allocateTranscriptCounts(transData, transDataList, startPos))
                        {
                            endOfTrans = true;
                            break;
                        }
                    }

                    if (endOfTrans)
                        break;
                }
            }

            // and generate fragments assuming an unspliced gene
            List<Integer> emptyTrans = Lists.newArrayList();

            if(commonExonicRegions.size() > 1)
            {
                long regionStart = geneReadData.getTranscriptsRange()[SE_START];
                long regionEnd = geneReadData.getTranscriptsRange()[SE_END];

                int exonicRegionIndex = 0;
                long currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];
                long nextExonicStart = commonExonicRegions.get(exonicRegionIndex + 1)[SE_START];

                for (long startPos = regionStart; startPos <= regionEnd - mCurrentFragSize; ++startPos)
                {
                    if (startPos <= currentExonicEnd)
                    {
                        // check possible transcript exonic matches
                        allocateUnsplicedCounts(transDataList, startPos);
                    }
                    else
                    {
                        // check for purely intronic fragments
                        if (startPos < nextExonicStart)
                        {
                            addTransComboData(UNSPLICED_ID, emptyTrans, UNSPLICED);
                        }
                        else
                        {
                            ++exonicRegionIndex;
                            currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];

                            if (exonicRegionIndex < commonExonicRegions.size() - 1)
                            {
                                nextExonicStart = commonExonicRegions.get(exonicRegionIndex + 1)[SE_START];
                            }
                            else
                            {
                                nextExonicStart = -1;
                            }
                        }
                    }
                }
            }
            else
            {
                // force an empty entry even though it won't have an category ratios set for it
                List<TranscriptComboData> emptyList = Lists.newArrayList(new TranscriptComboData(emptyTrans));
                mExpectedTransComboCounts.put(UNSPLICED_ID, emptyList);
            }
        }

        formTranscriptDefinitions(geneReadData);
    }

    private boolean allocateTranscriptCounts(final TranscriptData transData, final List<TranscriptData> transDataList, long startPos)
    {
        List<long[]> readRegions = Lists.newArrayList();
        List<long[]> spliceJunctions = Lists.newArrayList();

        FragmentMatchType matchType = generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        if(readRegions.isEmpty())
            return false;

        final List<Integer> longAndSplicedTrans = Lists.newArrayList();
        final List<Integer> shortTrans = Lists.newArrayList();

        if(matchType == SPLICED || matchType == LONG)
            longAndSplicedTrans.add(transData.TransId);
        else
            shortTrans.add(transData.TransId);

        // now check whether these regions are supported by each other transcript
        for(TranscriptData otherTransData : transDataList)
        {
            if(transData == otherTransData)
                continue;

            if(readsSupportFragment(otherTransData, readRegions, matchType, spliceJunctions))
            {
                if(matchType == SPLICED || matchType == LONG)
                    longAndSplicedTrans.add(otherTransData.TransId);
                else
                    shortTrans.add(otherTransData.TransId);
            }
        }

        if(!longAndSplicedTrans.isEmpty())
        {
            addTransComboData(transData.TransName, longAndSplicedTrans, SPLICED);
        }
        else
        {
            addTransComboData(transData.TransName, shortTrans, SHORT);
        }

        return true;
    }

    private void allocateUnsplicedCounts(final List<TranscriptData> transDataList, long startPos)
    {
        List<long[]> readRegions = Lists.newArrayList();
        List<long[]> noSpliceJunctions = Lists.newArrayList();

        // the unspliced case
        long firstReadEnd = startPos + mCurrentReadLength - 1;
        long secondReadEnd = startPos + mCurrentFragSize - 1;
        long secondReadStart = secondReadEnd - mCurrentReadLength + 1;

        if(firstReadEnd >= secondReadStart - 1)
        {
            // continuous reads so merge into one
            readRegions.add(new long[] {startPos, secondReadEnd});
        }
        else
        {
            readRegions.add(new long[] {startPos, firstReadEnd});
            readRegions.add(new long[] {secondReadStart, secondReadEnd});
        }

        final List<Integer> shortTrans = Lists.newArrayList();

        // check whether these unspliced reads support exonic regions
        for(TranscriptData transData : transDataList)
        {
            if(readsSupportFragment(transData, readRegions, SHORT, noSpliceJunctions))
            {
                shortTrans.add(transData.TransId);
            }
        }

        addTransComboData(UNSPLICED_ID, shortTrans, !shortTrans.isEmpty() ? SHORT : UNSPLICED);
    }

    private void addTransComboData(final String transId, final List<Integer> transcripts, FragmentMatchType transMatchType)
    {
        List<TranscriptComboData> transComboDataList = mExpectedTransComboCounts.get(transId);

        if(transComboDataList == null)
        {
            transComboDataList = Lists.newArrayList();
            mExpectedTransComboCounts.put(transId, transComboDataList);
        }

        TranscriptComboData matchingCounts = transComboDataList.stream()
                .filter(x -> x.matches(transcripts)).findFirst().orElse(null);

        if(matchingCounts == null)
        {
            matchingCounts = new TranscriptComboData(transcripts);
            transComboDataList.add(matchingCounts);
        }

        matchingCounts.addCounts(transMatchType, mCurrentFragFrequency);
    }

    public FragmentMatchType generateImpliedFragment(final TranscriptData transData, long startPos, List<long[]> readRegions, List<long[]> spliceJunctions)
    {
        readRegions.clear();
        spliceJunctions.clear();

        // set out the fragment reads either within a single exon or spanning one or more
        int exonCount = transData.exons().size();
        final ExonData lastExon = transData.exons().get(exonCount - 1);

        if(startPos + mCurrentFragSize - 1 > lastExon.ExonEnd)
            return UNSPLICED;

        FragmentMatchType matchType = SHORT;

        int remainingReadBases = mCurrentReadLength;
        boolean overlappingReads = (mCurrentFragSize - 2 * mCurrentReadLength) < 1;
        int remainingInterimBases = !overlappingReads ? mCurrentFragSize - 2 * mCurrentReadLength + 1 : 0;
        long nextRegionStart = startPos;
        int readsAdded = 0;

        for(int i = 0; i < exonCount; ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(nextRegionStart > exon.ExonEnd)
                continue;

            if(!readRegions.isEmpty())
            {
                if(matchType != SPLICED)
                    matchType = LONG;
            }

            if(readsAdded == 1 && remainingInterimBases > 0)
            {
                if(nextRegionStart + remainingInterimBases - 1 >= exon.ExonEnd)
                {
                    if(i >= exonCount - 1)
                    {
                        readRegions.clear();
                        return UNSPLICED;
                    }

                    nextRegionStart = transData.exons().get(i + 1).ExonStart;

                    remainingInterimBases -= exon.ExonEnd - exon.ExonStart + 1;
                    continue;

                }

                nextRegionStart += remainingInterimBases;
                remainingInterimBases = 0;
            }

            long regionEnd = min(nextRegionStart + remainingReadBases - 1, exon.ExonEnd);
            long regionLength = (int)(regionEnd - nextRegionStart + 1);
            remainingReadBases -= regionLength;
            readRegions.add(new long[] {nextRegionStart, regionEnd});
            boolean spansExonEnd = remainingReadBases > 0;

            if(remainingReadBases == 0)
            {
                ++readsAdded;

                if (readsAdded == 2)
                    break;

                if(overlappingReads)
                {
                    if(mCurrentFragSize <= mCurrentReadLength)
                        break;

                    remainingReadBases = mCurrentFragSize - mCurrentReadLength;
                }
                else
                {
                    remainingReadBases = mCurrentReadLength;
                }
            }

            // is the remainder of this exon long enough to match again?
            if(!overlappingReads)
                nextRegionStart = regionEnd + remainingInterimBases;
            else
                nextRegionStart = regionEnd + 1;

            if(regionEnd == exon.ExonEnd || nextRegionStart > exon.ExonEnd)
            {
                if(i == exonCount - 1)
                {
                    readRegions.clear();
                    return UNSPLICED;
                }

                // will move onto the next exon for further matching
                nextRegionStart = transData.exons().get(i + 1).ExonStart;

                if(spansExonEnd && regionEnd == exon.ExonEnd)
                {
                    matchType = SPLICED;
                    spliceJunctions.add(new long[] {exon.ExonEnd, nextRegionStart});
                }

                remainingInterimBases -= exon.ExonEnd - regionEnd;
                continue;
            }
            else
            {
                remainingInterimBases = 0;
            }

            // start the next match within this same exon
            regionEnd = min(nextRegionStart + remainingReadBases - 1, exon.ExonEnd);
            regionLength = (int)(regionEnd - nextRegionStart + 1);
            remainingReadBases -= regionLength;

            readRegions.add(new long[] { nextRegionStart, regionEnd });

            if(remainingReadBases > 0 && regionEnd == exon.ExonEnd)
                matchType = SPLICED;

            if(remainingReadBases == 0)
                break;

            if(i == exonCount - 1)
            {
                readRegions.clear();
                return UNSPLICED;
            }

            // will move onto the next exon for further matching
            nextRegionStart = transData.exons().get(i + 1).ExonStart;
        }

        // merge adjacent regions from overlapping reads
        if(overlappingReads)
        {
            int index = 0;
            while(index < readRegions.size() - 1)
            {
                long regionEnd = readRegions.get(index)[SE_END];
                long regionStart = readRegions.get(index + 1)[SE_START];

                if(regionStart == regionEnd + 1)
                {
                    readRegions.get(index)[SE_END] = readRegions.get(index + 1)[SE_END];
                    readRegions.remove(index + 1);
                    break;
                }

                ++index;
            }
        }

        return matchType;
    }

    public boolean readsSupportFragment(
            final TranscriptData transData, List<long[]> readRegions, FragmentMatchType requiredMatchType, List<long[]> spliceJunctions)
    {
        long regionsStart = readRegions.get(0)[SE_START];
        long regionsEnd = readRegions.get(readRegions.size() - 1)[SE_END];

        if(requiredMatchType == SHORT)
        {
            return (transData.exons().stream().anyMatch(x -> x.ExonStart <= regionsStart && x.ExonEnd >= regionsEnd));
        }
        else
        {
            // first check for matching splice junctions
            for(long[] spliceJunction : spliceJunctions)
            {
                long spliceStart = spliceJunction[SE_START];
                long spliceEnd = spliceJunction[SE_END];
                boolean matched = false;

                for (int i = 0; i < transData.exons().size() - 1; ++i)
                {
                    ExonData exon = transData.exons().get(i);
                    ExonData nextExon = transData.exons().get(i + 1);

                    if (exon.ExonEnd == spliceStart && nextExon.ExonStart == spliceEnd)
                    {
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    return false;
            }

            // none of the reads can breach an exon boundary or skip an exon
            int readIndex = 0;
            long regionStart = readRegions.get(readIndex)[SE_START];
            long regionEnd = readRegions.get(readIndex)[SE_END];
            int exonsMatched = 0;

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                ExonData exon = transData.exons().get(i);

                // region before the next exon even starts
                if(regionEnd < exon.ExonStart)
                    return false;

                // invalid if overlaps but not fully contained
                if(regionsStart >= exon.ExonStart && regionsStart <= exon.ExonEnd && regionEnd > exon.ExonEnd)
                    return false;
                else if(regionsEnd >= exon.ExonStart && regionsEnd <= exon.ExonEnd && regionStart < exon.ExonStart)
                    return false;

                ++exonsMatched;

                while(true)
                {
                    if (regionStart >= exon.ExonStart && regionEnd <= exon.ExonEnd)
                    {
                        ++readIndex;

                        if(readIndex >= readRegions.size())
                            break;

                        regionStart = readRegions.get(readIndex)[SE_START];
                        regionEnd = readRegions.get(readIndex)[SE_END];
                    }
                    else
                    {
                        // next region may match the next exon
                        break;
                    }
                }

                if(readIndex >= readRegions.size())
                    break;
            }

            return exonsMatched > 1;
        }
    }

    private void formTranscriptDefinitions(final GeneReadData geneReadData)
    {
        collectCategories();

        int categoryCount = mCurrentExpRatesData.Categories.size();

        mExpectedTransComboCounts.keySet().forEach(x -> mCurrentExpRatesData.TranscriptIds.add(x));

        mCurrentExpRatesData.initialiseTranscriptDefinitions();

        for(int transIndex = 0; transIndex < mCurrentExpRatesData.TranscriptIds.size(); ++transIndex)
        {
            final String transId = mCurrentExpRatesData.TranscriptIds.get(transIndex);

            double[] categoryCounts = new double[categoryCount];

            final List<TranscriptComboData> transCounts = mExpectedTransComboCounts.get(transId);

            for(TranscriptComboData tcData : transCounts)
            {
                final String transKey = !tcData.getTranscriptIds().isEmpty() ? tcData.getTranscriptsKey() : UNSPLICED_ID;

                final int[] counts = tcData.getCounts();

                for(int i = 0; i < counts.length; ++i)
                {
                    double catCount = counts[i];

                    if(catCount > 0)
                    {
                        final String categoryStr = formCategory(transKey, FragmentMatchType.intAsType(i));
                        int categoryId = mCurrentExpRatesData.getCategoryIndex(categoryStr);

                        if(categoryId < 0)
                        {
                            RE_LOGGER.error("invalid category index from transKey({})", transKey);
                            return;
                        }

                        // if(i != TC_SPLICED)
                        //    catCount *= mConfig.UnsplicedWeight;

                        categoryCounts[categoryId] = catCount;
                    }
                }
            }

            // convert counts to ratios and add against this transcript definition
            convertToPercentages(categoryCounts);
            mCurrentExpRatesData.getTranscriptDefinitions().setCol(transIndex, categoryCounts);
        }

        mCurrentExpRatesData.getTranscriptDefinitions().cacheTranspose();

        if(mConfig.GenerateExpectedRates)
            writeExpectedRates(mExpRateWriter, geneReadData,mCurrentExpRatesData);
    }

    private static final String UNSPLICED_STR = "UNSPLC";
    public static final int UNSPLICED_CAT_INDEX = 0;

    private void collectCategories()
    {
        mCurrentExpRatesData.addCategory(UNSPLICED_ID);

        for(Map.Entry<String,List<TranscriptComboData>> entry : mExpectedTransComboCounts.entrySet())
        {
            for(TranscriptComboData tcData : entry.getValue())
            {
                boolean hasTrans = !tcData.getTranscriptIds().isEmpty();
                final String transKey = hasTrans ? tcData.getTranscriptsKey() : UNSPLICED_ID;

                if(tcData.getShortCount() > 0)
                {
                    mCurrentExpRatesData.addCategory(formCategory(transKey, SHORT));
                }

                if(tcData.getSplicedCount() > 0)
                {
                    mCurrentExpRatesData.addCategory(transKey);
                }
            }
        }
    }

    public static String formCategory(final String transKey, FragmentMatchType countsType)
    {
        if(countsType == SHORT)
            return String.format("%s-%s", transKey, UNSPLICED_STR);
        else if(countsType == SPLICED)
            return transKey;
        else
            return UNSPLICED_ID;
    }

    public void setFragmentLengthData(int length, int frequency)
    {
        mCurrentFragSize = length;
        mCurrentFragFrequency = frequency;
        mCurrentReadLength = min(mConfig.ReadLength, mCurrentFragSize);
    }

    public static BufferedWriter createWriter(final RnaExpConfig config)
    {
        try
        {
            String outputFileName;

            if(config.ExpRatesFile != null)
            {
                outputFileName = config.ExpRatesFile;
            }
            else
            {
                outputFileName = config.OutputDir + EXP_RATES_FILE_ID;

                if(config.OutputIdentifier != null)
                    outputFileName += "." + config.OutputIdentifier;

                outputFileName += ".csv";
            }

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,TransId,Category,Rate");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            RE_LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeExpectedRates(
            final BufferedWriter writer, final GeneReadData geneReadData, final ExpectedRatesData expExpData)
    {
        if(writer == null)
            return;

        if(expExpData.getTranscriptDefinitions() == null || expExpData.Categories.isEmpty() || expExpData.TranscriptIds.isEmpty())
            return;

        try
        {
            final SigMatrix transcriptDefinitions = expExpData.getTranscriptDefinitions();
            final List<String> categories = expExpData.Categories;
            final List<String> transcriptIds = expExpData.TranscriptIds;

            final String geneId = geneReadData.GeneData.GeneId;
            final String geneName = geneReadData.GeneData.GeneName;

            for(int i = 0; i < transcriptDefinitions.Cols; ++i)
            {
                final String transcriptId = transcriptIds.get(i);

                for(int j = 0; j < transcriptDefinitions.Rows; ++j)
                {
                    final String category = categories.get(j);

                    double expRate = transcriptDefinitions.get(j, i);

                    if(expRate < 0.0001)
                        continue;

                    writer.write(String.format("%s,%s,%s,%s,%.4f",
                            geneId, geneName, transcriptId, category, expRate));
                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
        }
    }


}
