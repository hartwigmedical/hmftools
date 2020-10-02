package com.hartwig.hmftools.isofox.expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.LONG;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.UNSPLICED;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_TRANS_COUNTS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.common.sigs.SigMatrix;

public class ExpectedRatesGenerator
{
    private final IsofoxConfig mConfig;

    // map of transcript (or unspliced) to all expected category counts covering it (and others)
    private final Map<String,List<CategoryCountsData>> mTransCategoryCounts;

    private GeneCollection mGeneCollection;
    private ExpectedRatesData mCurrentExpRatesData;

    private int mCurrentFragSize;
    private int mFragSizeIndex;
    private int mCurrentFragFrequency;
    private int mReadLength;

    private final List<ExpectedRatesData> mExpectedRatesDataList;

    private final BufferedWriter mExpRateWriter;

    public static final int FL_LENGTH = 0;
    public static final int FL_FREQUENCY = 1;

    public ExpectedRatesGenerator(final IsofoxConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mCurrentFragSize = 0;
        mFragSizeIndex = 0;
        mCurrentFragFrequency = 0;
        mReadLength = mConfig.ReadLength;

        mTransCategoryCounts = Maps.newHashMap();
        mCurrentExpRatesData = null;
        mGeneCollection = null;
        mExpectedRatesDataList = Lists.newArrayList();

        mExpRateWriter = resultsWriter != null ? resultsWriter.getExpRatesWriter() : null;
    }

    public static ExpectedRatesGenerator from(final IsofoxConfig config)
    {
        return new ExpectedRatesGenerator(config, null);
    }

    public Map<String,List<CategoryCountsData>> getTransComboData() { return mTransCategoryCounts; }
    public ExpectedRatesData getExpectedRatesData() { return mCurrentExpRatesData; }

    public void generateExpectedRates(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;
        mTransCategoryCounts.clear();
        mCurrentExpRatesData = new ExpectedRatesData(mGeneCollection.chrId());

        final List<int[]> commonExonicRegions = mGeneCollection.getCommonExonicRegions();

        // apply fragment reads across each transcript as though it were fully transcribed
        final List<TranscriptData> transDataList = mGeneCollection.getTranscripts();

        for(mFragSizeIndex = 0; mFragSizeIndex < mConfig.FragmentSizeData.size(); ++mFragSizeIndex)
        {
            final FragmentSize flData = mConfig.FragmentSizeData.get(mFragSizeIndex);
            mCurrentFragSize = flData.Length;
            mCurrentFragFrequency = flData.Frequency;

            for (TranscriptData transData : transDataList)
            {
                boolean endOfTrans = false;
                List<TranscriptData> candidateTrans = Lists.newArrayList(transDataList);
                candidateTrans.remove(transData);

                for (ExonData exon : transData.exons())
                {
                    for (int startPos = exon.ExonStart; startPos <= exon.ExonEnd; ++startPos)
                    {
                        cullTranscripts(candidateTrans, startPos);

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
                int regionStart = mGeneCollection.regionBounds()[SE_START];
                int regionEnd = mGeneCollection.regionBounds()[SE_END];

                int exonicRegionIndex = 0;
                int currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];
                int nextExonicStart = commonExonicRegions.get(exonicRegionIndex + 1)[SE_START];

                List<TranscriptData> candidateTrans = Lists.newArrayList(transDataList);

                for (int startPos = regionStart; startPos <= regionEnd - mCurrentFragSize; ++startPos)
                {
                    final List<String> unsplicedGenes = findUnsplicedGenes(startPos);

                    // cull the set of possible transcripts
                    cullTranscripts(candidateTrans, startPos);

                    if (startPos <= currentExonicEnd)
                    {
                        // check possible transcript exonic matches
                        allocateUnsplicedCounts(transDataList, startPos, unsplicedGenes);
                    }
                    else
                    {
                        // check for purely intronic fragments
                        if (startPos < nextExonicStart)
                        {
                            addUnsplicedCountsData(emptyTrans, unsplicedGenes);
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
                // force an empty entry even though it won't have any category ratios set for it
                List<String> allGeneIds = mGeneCollection.genes().stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());
                CategoryCountsData genesWithoutCounts = new CategoryCountsData(emptyTrans, allGeneIds);
                genesWithoutCounts.initialiseLengthCounts(mConfig.FragmentSizeData.size());
                List<CategoryCountsData> emptyList = Lists.newArrayList(genesWithoutCounts);

                for(GeneReadData gene : mGeneCollection.genes())
                {
                    mTransCategoryCounts.put(gene.GeneData.GeneId, emptyList);
                }
            }
        }

        // add in any genes which ended up without counts, ie those with a single exon
        for(GeneReadData gene : geneCollection.genes())
        {
            final String geneId = gene.GeneData.GeneId;
            if(mTransCategoryCounts.containsKey(geneId))
                continue;

            CategoryCountsData genesWithoutCounts = new CategoryCountsData(Lists.newArrayList(), Lists.newArrayList(geneId));
            genesWithoutCounts.initialiseLengthCounts(mConfig.FragmentSizeData.size());
            List<CategoryCountsData> emptyList = Lists.newArrayList(genesWithoutCounts);

            mTransCategoryCounts.put(geneId, emptyList);
        }

        if(mConfig.runFunction(EXPECTED_TRANS_COUNTS))
        {
            writeExpectedCounts(mExpRateWriter, geneCollection.chrId(), mTransCategoryCounts);
        }
        else
        {
            formTranscriptDefinitions(mTransCategoryCounts, mCurrentExpRatesData);

            if(mConfig.ApplyExpectedRates)
                mCurrentExpRatesData.getTranscriptDefinitions().cacheTranspose();

            if(mConfig.WriteExpectedRates)
            {
                writeExpectedRates(mExpRateWriter, mCurrentExpRatesData);
                // cache and write later
                // mExpectedRatesDataList.add(mCurrentExpRatesData);
            }
        }
    }

    private void cullTranscripts(final List<TranscriptData> transcripts, int startPos)
    {
        int index = 0;
        while(index < transcripts.size())
        {
            if(transcripts.get(index).TransEnd < startPos)
                transcripts.remove(index);
            else
                ++index;
        }
    }

    private List<String> findUnsplicedGenes(int fragStart)
    {
        if(mGeneCollection.genes().size() == 1)
        {
            if(mGeneCollection.getExonRegions().size() == 1)
                return Lists.newArrayList();
            else
                return Lists.newArrayList(mGeneCollection.genes().get(0).GeneData.GeneId);
        }

        int fragEnd = fragStart + mCurrentFragSize - 1;

        return mGeneCollection.genes().stream()
                .filter(x -> x.getExonRegions().size() > 1) // must have at least one unspliced region
                .filter(x -> positionsWithin(fragStart,fragEnd, x.GeneData.GeneStart, x.GeneData.GeneEnd))
                .map(x -> x.GeneData.GeneId).collect(Collectors.toList());
    }

    private boolean allocateTranscriptCounts(final TranscriptData transData, final List<TranscriptData> transDataList, int startPos)
    {
        List<int[]> readRegions = Lists.newArrayList();
        List<int[]> spliceJunctions = Lists.newArrayList();

        FragmentMatchType matchType = generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        if(readRegions.isEmpty())
            return false;

        final List<Integer> longAndSplicedTrans = Lists.newArrayList();
        final List<Integer> shortTrans = Lists.newArrayList();

        if(matchType == SPLICED || matchType == LONG)
            longAndSplicedTrans.add(transData.TransId);
        else
            shortTrans.add(transData.TransId);

        // now check whether these regions are supported by each other's transcript
        for(TranscriptData otherTransData : transDataList)
        {
            if(transData == otherTransData)
                continue;

            if(readsSupportTranscript(otherTransData, readRegions, matchType, spliceJunctions))
            {
                if(matchType == SPLICED || matchType == LONG)
                    longAndSplicedTrans.add(otherTransData.TransId);
                else
                    shortTrans.add(otherTransData.TransId);
            }
        }

        if(!longAndSplicedTrans.isEmpty())
        {
            addCountsData(transData.TransName, longAndSplicedTrans, Lists.newArrayList());
        }
        else
        {
            List<String> unsplicedGenes = findUnsplicedGenes(startPos);
            addCountsData(transData.TransName, shortTrans, unsplicedGenes);
        }

        return true;
    }

    private void allocateUnsplicedCounts(final List<TranscriptData> transDataList, int startPos, final List<String> unsplicedGenes)
    {
        List<int[]> readRegions = Lists.newArrayList();
        List<int[]> noSpliceJunctions = Lists.newArrayList();

        // the unspliced case
        int firstReadEnd = startPos + mReadLength - 1;
        int secondReadEnd = startPos + mCurrentFragSize - 1;
        int secondReadStart = secondReadEnd - mReadLength + 1;

        if(firstReadEnd >= secondReadStart - 1)
        {
            // continuous reads so merge into one
            readRegions.add(new int[] {startPos, secondReadEnd});
        }
        else
        {
            readRegions.add(new int[] {startPos, firstReadEnd});
            readRegions.add(new int[] {secondReadStart, secondReadEnd});
        }

        final List<Integer> shortTrans = Lists.newArrayList();

        // check whether these unspliced reads support exonic regions
        for(TranscriptData transData : transDataList)
        {
            if(readsSupportTranscript(transData, readRegions, SHORT, noSpliceJunctions))
            {
                shortTrans.add(transData.TransId);
            }
        }

        addUnsplicedCountsData(shortTrans, unsplicedGenes);
    }

    private void addUnsplicedCountsData(final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        unsplicedGenes.forEach(x -> addCountsData(x, transcripts, unsplicedGenes));
    }

    private void addCountsData(
            final String transName, final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        List<CategoryCountsData> transComboDataList = mTransCategoryCounts.get(transName);

        if(transComboDataList == null)
        {
            transComboDataList = Lists.newArrayList();
            mTransCategoryCounts.put(transName, transComboDataList);
        }

        CategoryCountsData matchingCounts = transComboDataList.stream()
                .filter(x -> x.matches(transcripts, unsplicedGenes)).findFirst().orElse(null);

        if(matchingCounts == null)
        {
            matchingCounts = new CategoryCountsData(transcripts, unsplicedGenes);

            if(mConfig.runFunction(EXPECTED_TRANS_COUNTS))
                matchingCounts.initialiseLengthCounts(mConfig.FragmentSizeData.size());

            transComboDataList.add(matchingCounts);
        }

        if(mConfig.runFunction(EXPECTED_TRANS_COUNTS))
            matchingCounts.addFragLengthCounts(mCurrentFragFrequency, mFragSizeIndex);
        else
            matchingCounts.addCounts(mCurrentFragFrequency);
    }

    public FragmentMatchType generateImpliedFragment(
            final TranscriptData transData, int startPos, List<int[]> readRegions, List<int[]> spliceJunctions)
    {
        readRegions.clear();
        spliceJunctions.clear();

        // set out the fragment reads either within a single exon or spanning one or more
        int exonCount = transData.exons().size();
        final ExonData lastExon = transData.exons().get(exonCount - 1);

        if(startPos + mCurrentFragSize - 1 > lastExon.ExonEnd)
            return UNSPLICED;

        FragmentMatchType matchType = SHORT;

        int remainingReadBases = mReadLength;
        boolean overlappingReads = (mCurrentFragSize - 2 * mReadLength) < 1;
        int remainingInterimBases = !overlappingReads ? mCurrentFragSize - 2 * mReadLength + 1 : 0;
        int nextRegionStart = startPos;
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

            int regionEnd = min(nextRegionStart + remainingReadBases - 1, exon.ExonEnd);
            int regionLength = (int)(regionEnd - nextRegionStart + 1);
            remainingReadBases -= regionLength;
            readRegions.add(new int[] {nextRegionStart, regionEnd});
            boolean spansExonEnd = remainingReadBases > 0;

            if(remainingReadBases == 0)
            {
                ++readsAdded;

                if (readsAdded == 2)
                    break;

                if(overlappingReads)
                {
                    if(mCurrentFragSize <= mReadLength)
                        break;

                    remainingReadBases = mCurrentFragSize - mReadLength;
                }
                else
                {
                    remainingReadBases = mReadLength;
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
                    spliceJunctions.add(new int[] {exon.ExonEnd, nextRegionStart});
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

            readRegions.add(new int[] { nextRegionStart, regionEnd });

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
                int regionEnd = readRegions.get(index)[SE_END];
                int regionStart = readRegions.get(index + 1)[SE_START];

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

    public boolean readsSupportTranscript(
            final TranscriptData transData, List<int[]> readRegions, FragmentMatchType requiredMatchType, List<int[]> spliceJunctions)
    {
        int regionsStart = readRegions.get(0)[SE_START];
        int regionsEnd = readRegions.get(readRegions.size() - 1)[SE_END];

        if(!positionsOverlap(regionsStart, regionsEnd, transData.TransStart, transData.TransEnd))
            return false;

        if(requiredMatchType == SHORT)
        {
            return (transData.exons().stream().anyMatch(x -> x.ExonStart <= regionsStart && x.ExonEnd >= regionsEnd));
        }
        else
        {
            // first check for matching splice junctions
            for(int[] spliceJunction : spliceJunctions)
            {
                int spliceStart = spliceJunction[SE_START];
                int spliceEnd = spliceJunction[SE_END];
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
            int regionStart = readRegions.get(readIndex)[SE_START];
            int regionEnd = readRegions.get(readIndex)[SE_END];
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

    public static void formTranscriptDefinitions(final Map<String,List<CategoryCountsData>> transCategoryCounts, ExpectedRatesData expRatesData)
    {
        // convert fragment counts in each category per transcript into the equivalent of a signature per transcript
        collectCategories(transCategoryCounts, expRatesData);

        int categoryCount = expRatesData.Categories.size();

        transCategoryCounts.keySet().forEach(x -> expRatesData.TranscriptIds.add(x));

        expRatesData.initialiseTranscriptDefinitions();

        for(int transIndex = 0; transIndex < expRatesData.TranscriptIds.size(); ++transIndex)
        {
            final String transId = expRatesData.TranscriptIds.get(transIndex);

            double[] categoryCounts = new double[categoryCount];

            final List<CategoryCountsData> transCounts = transCategoryCounts.get(transId);

            for(CategoryCountsData tcData : transCounts)
            {
                final String transKey = tcData.combinedKey();
                double fragmentCount = tcData.fragmentCount();

                if(fragmentCount > 0)
                {
                    int categoryId = expRatesData.getCategoryIndex(transKey);

                    if(categoryId < 0)
                    {
                        ISF_LOGGER.error("invalid category index from transKey({})", transKey);
                        return;
                    }

                    categoryCounts[categoryId] = fragmentCount;
                }
            }

            // convert counts to ratios and add against this transcript definition
            convertToPercentages(categoryCounts);
            expRatesData.getTranscriptDefinitions().setCol(transIndex, categoryCounts);
        }
    }

    public static void collectCategories(final Map<String,List<CategoryCountsData>> transCategoryCounts, ExpectedRatesData expRatesData)
    {
        for(Map.Entry<String,List<CategoryCountsData>> entry : transCategoryCounts.entrySet())
        {
            for(CategoryCountsData tcData : entry.getValue())
            {
                final String transKey = tcData.combinedKey();

                if(tcData.fragmentCount() > 0 || tcData.transcriptIds().isEmpty()) // force inclusion of unspliced gene categories
                {
                    expRatesData.addCategory(transKey);
                }
            }
        }
    }

    public void setFragmentLengthData(int length, int frequency)
    {
        mCurrentFragSize = length;
        mCurrentFragFrequency = frequency;
        mReadLength = min(mConfig.ReadLength, mCurrentFragSize);
    }

    public static final String EXP_COUNT_LENGTH_HEADER = "Length_";

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            String outputFileName;

            if(config.runFunction(EXPECTED_TRANS_COUNTS))
            {
                outputFileName = String.format("%sread_%d_%s", config.OutputDir, config.ReadLength, "exp_counts.csv");
            }
            else
            {
                outputFileName = config.formOutputFile("exp_rates.csv");
            }

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(config.runFunction(EXPECTED_TRANS_COUNTS))
            {
                writer.write("GeneSetId,TransId,Category");

                for(FragmentSize fragLength : config.FragmentSizeData)
                {
                    writer.write(String.format(",%s_%d", EXP_COUNT_LENGTH_HEADER, fragLength.Length));
                }
            }
            else
            {
                writer.write("GeneSetId,TransId,Category,Rate");
            }

            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
            return null;
        }
    }

    public void writeExpectedRatesData()
    {
        mExpectedRatesDataList.forEach(x -> writeExpectedRates(mExpRateWriter, x));
    }

    private synchronized static void writeExpectedCounts(
            final BufferedWriter writer, final String collectionId, final Map<String,List<CategoryCountsData>> transCountsMap)
    {
        if(writer == null)
            return;

        try
        {
            for(Map.Entry<String,List<CategoryCountsData>> entry : transCountsMap.entrySet())
            {
                final String transId = entry.getKey();
                final List<CategoryCountsData> countsList = entry.getValue();

                for(CategoryCountsData tcData : countsList)
                {
                    final int[] lengthCounts = tcData.fragmentCountsByLength();

                    writer.write(String.format("%s,%s,%s,%d", collectionId, transId, tcData.combinedKey(), lengthCounts[0]));

                    for (int i = 1; i < lengthCounts.length; ++i)
                    {
                        writer.write(String.format(",%d", lengthCounts[i]));
                    }

                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected counts file: {}", e.toString());
        }

    }

    public synchronized static void writeExpectedRates(
            final BufferedWriter writer, final ExpectedRatesData expExpData)
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

            for(int i = 0; i < transcriptDefinitions.Cols; ++i)
            {
                final String transcriptId = transcriptIds.get(i);

                for(int j = 0; j < transcriptDefinitions.Rows; ++j)
                {
                    final String category = categories.get(j);

                    double expRate = transcriptDefinitions.get(j, i);

                    if(expRate < 0.0001)
                        continue;

                    writer.write(String.format("%s,%s,%s,%.4f",
                            expExpData.Id, transcriptId, category, expRate));
                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
        }
    }


}
