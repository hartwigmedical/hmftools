package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.convertToPercentages;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.LONG;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.UNSPLICED;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExpectedExpressionRates
{
    private final RnaExpConfig mConfig;

    private final Map<String,List<TranscriptComboData>> mTransComboData;

    private final List<String> mCategories; // equivalent of buckets - 0-N transcripts and the fragment type (eg SHORT, SPLICED etc)
    private final List<String> mTranscriptIds; // equivalent of signature names - all transcripts and an UNSPLICED definition
    private SigMatrix mTranscriptDefinitions;

    private int mCurrentFragSize;
    private int mCurrentFragFrequency;
    private int mCurrentReadLength;

    private BufferedWriter mExpRateWriter;

    public static final String UNSPLICED_ID = "UNSPLICED";

    public static final int FL_SIZE = 0;
    public static final int FL_FREQUENCY = 1;

    private static final Logger LOGGER = LogManager.getLogger(ExpectedExpressionRates.class);

    public ExpectedExpressionRates(final RnaExpConfig config)
    {
        mConfig = config;
        mCurrentFragSize = 0;
        mCurrentFragFrequency = 0;
        mCurrentReadLength = mConfig.ReadLength;

        mTransComboData = Maps.newHashMap();
        mCategories = Lists.newArrayList();
        mTranscriptIds = Lists.newArrayList();
        mTranscriptDefinitions = null;
        mExpRateWriter = null;
    }

    public Map<String,List<TranscriptComboData>> getTransComboData() { return mTransComboData; }

    public List<String> getCategories() { return mCategories; }
    public List<String> getTranscriptNames() { return mTranscriptIds; }
    public SigMatrix getTranscriptDefinitions() { return mTranscriptDefinitions; }

    public void setFragmentLengthData(int length, int frequency)
    {
        mCurrentFragSize = length;
        mCurrentFragFrequency = frequency;
        mCurrentReadLength = min(mConfig.ReadLength, mCurrentFragSize);
    }

    public boolean validData()
    {
        if(mCategories.isEmpty() || mTranscriptDefinitions == null || mTransComboData.isEmpty())
            return false;

        if(mTranscriptDefinitions.Cols != mTranscriptIds.size())
            return false;

        if(mTranscriptDefinitions.Rows != mCategories.size())
            return false;

        return true;
    }

    public void generateExpectedRates(final GeneReadData geneReadData)
    {
        clearState();

        final List<long[]> commonExonicRegions = geneReadData.getCommonExonicRegions();

        if(commonExonicRegions.size() < 2)
            return;

        // process each transcript as though it were transcribed
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
            long regionStart = geneReadData.getTranscriptsRange()[SE_START];
            long regionEnd = geneReadData.getTranscriptsRange()[SE_END];

            int exonicRegionIndex = 0;
            long currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];
            long nextExonicStart = commonExonicRegions.get(exonicRegionIndex + 1)[SE_START];

            List<String> emptyTrans = Lists.newArrayList();

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

        formTranscriptDefinitions(geneReadData);
    }

    private boolean allocateTranscriptCounts(final TranscriptData transData, final List<TranscriptData> transDataList, long startPos)
    {
        List<long[]> readRegions = Lists.newArrayList();
        List<long[]> spliceJunctions = Lists.newArrayList();

        FragmentMatchType matchType = generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        if(readRegions.isEmpty())
            return false;

        final List<String> longAndSplicedTrans = Lists.newArrayList();
        final List<String> shortTrans = Lists.newArrayList();

        if(matchType == SPLICED || matchType == LONG)
            longAndSplicedTrans.add(transData.TransName);
        else
            shortTrans.add(transData.TransName);

        // now check whether these regions are supported by each other transcript
        for(TranscriptData otherTransData : transDataList)
        {
            if(transData == otherTransData)
                continue;

            if(readsSupportFragment(otherTransData, readRegions, matchType, spliceJunctions))
            {
                if(matchType == SPLICED || matchType == LONG)
                    longAndSplicedTrans.add(otherTransData.TransName);
                else
                    shortTrans.add(otherTransData.TransName);
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

        final List<String> shortTrans = Lists.newArrayList();

        // check whether these unspliced reads support exonic regions
        for(TranscriptData transData : transDataList)
        {
            if(readsSupportFragment(transData, readRegions, SHORT, noSpliceJunctions))
            {
                shortTrans.add(transData.TransName);
            }
        }

        addTransComboData(UNSPLICED_ID, shortTrans, !shortTrans.isEmpty() ? SHORT : UNSPLICED);
    }

    private void addTransComboData(final String transId, final List<String> transcripts, FragmentMatchType transMatchType)
    {
        List<TranscriptComboData> transComboDataList = mTransComboData.get(transId);

        if(transComboDataList == null)
        {
            transComboDataList = Lists.newArrayList();
            mTransComboData.put(transId, transComboDataList);
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

        int categoryCount = mCategories.size();
        int definitionsCount = mTransComboData.size();

        mTranscriptDefinitions = new SigMatrix(categoryCount, definitionsCount);

        for(Map.Entry<String,List<TranscriptComboData>> entry : mTransComboData.entrySet())
        {
            final String transId = entry.getKey();

            int definitionIndex = mTranscriptIds.size();
            mTranscriptIds.add(transId);

            double[] categoryCounts = new double[categoryCount];

            for(TranscriptComboData tcData : entry.getValue())
            {
                final String transKey = !tcData.getTranscripts().isEmpty() ? tcData.getTranscriptsKey() : UNSPLICED_ID;

                final int[] counts = tcData.getCounts();

                for(int i = 0; i < counts.length; ++i)
                {
                    double catCount = counts[i];

                    if(catCount > 0)
                    {
                        final String categoryStr = formCategory(transKey, FragmentMatchType.intAsType(i));
                        int categoryId = getCategoryIndex(categoryStr);

                        if(categoryId < 0)
                        {
                            LOGGER.error("invalid category index from transKey({})", transKey);
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
            mTranscriptDefinitions.setCol(definitionIndex, categoryCounts);
        }

        mTranscriptDefinitions.cacheTranspose();

        if(mConfig.WriteExpectedRates)
            writeExpectedRates(geneReadData);
    }

    public double[] generateTranscriptCounts(final GeneReadData geneReadData, final List<TranscriptComboData> transComboData)
    {
        double[] categoryCounts = new double[mCategories.size()];

        int skippedComboCounts = 0;

        for(TranscriptComboData tcData : transComboData)
        {
            final String transKey = !tcData.getTranscripts().isEmpty() ? tcData.getTranscriptsKey() : UNSPLICED_ID;

            int shortCount = tcData.getShortCount();
            int splicedCount = tcData.getSplicedCount();

            if(shortCount > 0)
            {
                final String categoryStr = formCategory(transKey, SHORT);
                int categoryId = getCategoryIndex(categoryStr);

                // for now if a category isn't found just log and then ignore the count in it
                if(categoryId < 0)
                {
                    LOGGER.debug("category({}) skipped with count({})", categoryStr, shortCount);
                    skippedComboCounts += shortCount;
                }
                else
                {
                    categoryCounts[categoryId] = shortCount; //  * mConfig.UnsplicedWeight;
                }
            }

            if(splicedCount > 0)
            {
                final String categoryStr = formCategory(transKey, SPLICED);
                int categoryId = getCategoryIndex(categoryStr);

                if(categoryId < 0)
                {
                    LOGGER.debug("category({}) skipped with count({})", categoryStr, splicedCount);
                    skippedComboCounts += splicedCount;
                }
                else
                {
                    categoryCounts[categoryId] = splicedCount;
                }
            }
        }

        if(skippedComboCounts > 0)
        {
            double totalCounts = sumVector(categoryCounts) + skippedComboCounts;

            LOGGER.info(String.format("gene(%s) skippedCounts(%d perc=%.3f of total=%.0f)",
                    geneReadData.GeneData.GeneName, skippedComboCounts, skippedComboCounts/totalCounts, totalCounts));
        }

        return categoryCounts;
    }

    private int getCategoryIndex(final String category)
    {
        for(int i = 0; i < mCategories.size(); ++i)
        {
            if(mCategories.get(i).equals(category))
                return i;
        }

        return -1;
    }

    private void addCategory(final String category)
    {
        if(!mCategories.contains(category))
            mCategories.add(category);
    }

    private static final String UNSPLICED_STR = "UNSPLC";
    public static final int UNSPLICED_CAT_INDEX = 0;

    private void collectCategories()
    {
        addCategory(UNSPLICED_ID);

        for(Map.Entry<String,List<TranscriptComboData>> entry : mTransComboData.entrySet())
        {
            for(TranscriptComboData tcData : entry.getValue())
            {
                boolean hasTrans = !tcData.getTranscripts().isEmpty();
                final String transKey = hasTrans ? tcData.getTranscriptsKey() : UNSPLICED_ID;

                if(tcData.getShortCount() > 0)
                {
                    addCategory(formCategory(transKey, SHORT));
                }

                if(tcData.getSplicedCount() > 0)
                {
                    addCategory(transKey);
                }
            }
        }
    }

    private static String formCategory(final String transKey, FragmentMatchType countsType)
    {
        if(countsType == SHORT)
            return String.format("%s-%s", transKey, UNSPLICED_STR);
        else if(countsType == SPLICED)
            return transKey;
        else
            return UNSPLICED_ID;
    }

    private void clearState()
    {
        mTransComboData.clear();
        mTranscriptDefinitions = null;
        mCategories.clear();
        mTranscriptIds.clear();
    }

    private void writeExpectedRates(final GeneReadData geneReadData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        if(mTranscriptDefinitions == null || mCategories.isEmpty() || mTranscriptIds.isEmpty())
            return;

        try
        {
            if(mExpRateWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_EXP_TRAN_RATES.csv";

                mExpRateWriter = createBufferedWriter(outputFileName, false);
                mExpRateWriter.write("GeneId,GeneName,Transcript,Category,Rate");
                mExpRateWriter.newLine();
            }

            final String geneId = geneReadData.GeneData.GeneId;
            final String geneName = geneReadData.GeneData.GeneName;

            for(int i = 0; i < mTranscriptDefinitions.Cols; ++i)
            {
                final String transcriptId = mTranscriptIds.get(i);

                for(int j = 0; j < mTranscriptDefinitions.Rows; ++j)
                {
                    final String category = mCategories.get(j);

                    double expRate = mTranscriptDefinitions.get(j, i);

                    if(expRate < 0.0001)
                        continue;

                    mExpRateWriter.write(String.format("%s,%s,%s,%s,%.4f",
                            geneId, geneName, transcriptId, category, expRate));
                    mExpRateWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mExpRateWriter);
    }


}
