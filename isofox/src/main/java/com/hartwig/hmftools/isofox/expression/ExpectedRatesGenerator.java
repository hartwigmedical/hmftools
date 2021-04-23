package com.hartwig.hmftools.isofox.expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.BamFragmentAllocator.calcFragmentLength;
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
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.common.utils.Matrix;

public class ExpectedRatesGenerator
{
    private final IsofoxConfig mConfig;

    // map of transcript or unspliced gene to all expected category counts covering it (and others)
    private final Map<String,List<CategoryCountsData>> mTransCategoryCountsMap;

    // consolidated list of expected category counts
    private final List<CategoryCountsData> mTransCategoryCounts;

    private GeneCollection mGeneCollection;
    private ExpectedRatesData mCurrentExpRatesData;

    private int mCurrentFragSize;
    private int mFragSizeIndex;
    private int mCurrentFragFrequency;
    private int mReadLength;

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

        mTransCategoryCountsMap = Maps.newHashMap();
        mTransCategoryCounts = Lists.newArrayList();
        mCurrentExpRatesData = null;
        mGeneCollection = null;

        mExpRateWriter = resultsWriter != null ? resultsWriter.getExpRatesWriter() : null;
    }

    public static ExpectedRatesGenerator from(final IsofoxConfig config)
    {
        return new ExpectedRatesGenerator(config, null);
    }

    public List<CategoryCountsData> getTransComboData() { return mTransCategoryCounts; }
    public ExpectedRatesData getExpectedRatesData() { return mCurrentExpRatesData; }

    public void generateExpectedRates(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;
        mTransCategoryCountsMap.clear();
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
                    for (int startPos = exon.Start; startPos <= exon.End; ++startPos)
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
                    mTransCategoryCountsMap.put(gene.GeneData.GeneId, emptyList);
                }
            }
        }

        // add in any genes which ended up without counts, ie those with a single exon
        for(GeneReadData gene : geneCollection.genes())
        {
            final String geneId = gene.GeneData.GeneId;

            if(mTransCategoryCountsMap.containsKey(geneId))
                continue;

            CategoryCountsData genesWithoutCounts = new CategoryCountsData(Lists.newArrayList(), Lists.newArrayList(geneId));
            genesWithoutCounts.initialiseLengthCounts(mConfig.FragmentSizeData.size());
            List<CategoryCountsData> emptyList = Lists.newArrayList(genesWithoutCounts);

            mTransCategoryCountsMap.put(geneId, emptyList);
        }

        buildUniqueCategoryCounts();

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
            }
        }
    }

    private void buildUniqueCategoryCounts()
    {
        // take the category data across all transcripts and genes and convert it into a unique list
        for(List<CategoryCountsData> categoryCounts : mTransCategoryCountsMap.values())
        {
            for(CategoryCountsData catCountsData : categoryCounts)
            {
                CategoryCountsData matchedData = mTransCategoryCounts.stream()
                        .filter(x -> x.combinedKey().equals(catCountsData.combinedKey())
                                || x.matches(catCountsData.transcriptIds(), catCountsData.unsplicedGeneIds()))
                        .findFirst().orElse(null);

                if(matchedData == null)
                {
                    mTransCategoryCounts.add(catCountsData);
                }
                else
                {
                    matchedData.mergeCounts(catCountsData);
                }
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
            if(mGeneCollection.genes().get(0).hasUnsplicedRegions())
                return Lists.newArrayList(mGeneCollection.genes().get(0).GeneData.GeneId);
            else
                return Lists.newArrayList();
        }

        int fragEnd = fragStart + mCurrentFragSize - 1;

        return mGeneCollection.genes().stream()
                .filter(x -> x.hasUnsplicedRegions()) // must have at least one unspliced region
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

    private void addCountsData(final String transName, final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        List<CategoryCountsData> transComboDataList = mTransCategoryCountsMap.get(transName);

        if(transComboDataList == null)
        {
            transComboDataList = Lists.newArrayList();
            mTransCategoryCountsMap.put(transName, transComboDataList);
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

        if(startPos + mCurrentFragSize - 1 > lastExon.End)
            return UNSPLICED;

        FragmentMatchType matchType = SHORT;

        int remainingReadBases = min(mReadLength, mCurrentFragSize);
        boolean overlappingReads = (mCurrentFragSize - 2 * mReadLength) < 1;
        int remainingInterimBases = !overlappingReads ? mCurrentFragSize - 2 * mReadLength + 1 : 0;
        int nextRegionStart = startPos;
        int readsAdded = 0;

        for(int i = 0; i < exonCount; ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(nextRegionStart > exon.End)
                continue;

            if(!readRegions.isEmpty())
            {
                if(matchType != SPLICED)
                    matchType = LONG;
            }

            if(readsAdded == 1 && remainingInterimBases > 0)
            {
                if(nextRegionStart + remainingInterimBases - 1 >= exon.End)
                {
                    if(i >= exonCount - 1)
                    {
                        readRegions.clear();
                        return UNSPLICED;
                    }

                    nextRegionStart = transData.exons().get(i + 1).Start;

                    remainingInterimBases -= exon.End - exon.Start + 1;
                    continue;

                }

                nextRegionStart += remainingInterimBases;
                remainingInterimBases = 0;
            }

            int regionEnd = min(nextRegionStart + remainingReadBases - 1, exon.End);
            int regionLength = regionEnd - nextRegionStart + 1;
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

            if(regionEnd == exon.End || nextRegionStart > exon.End)
            {
                if(i == exonCount - 1)
                {
                    readRegions.clear();
                    return UNSPLICED;
                }

                // will move onto the next exon for further matching
                nextRegionStart = transData.exons().get(i + 1).Start;

                if(spansExonEnd && regionEnd == exon.End)
                {
                    matchType = SPLICED;
                    spliceJunctions.add(new int[] {exon.End, nextRegionStart});
                }

                remainingInterimBases -= exon.End - regionEnd;
                continue;
            }
            else
            {
                remainingInterimBases = 0;
            }

            // start the next match within this same exon
            regionEnd = min(nextRegionStart + remainingReadBases - 1, exon.End);
            regionLength = (int)(regionEnd - nextRegionStart + 1);
            remainingReadBases -= regionLength;

            readRegions.add(new int[] { nextRegionStart, regionEnd });

            if(remainingReadBases > 0 && regionEnd == exon.End)
                matchType = SPLICED;

            if(remainingReadBases == 0)
                break;

            if(i == exonCount - 1)
            {
                readRegions.clear();
                return UNSPLICED;
            }

            // will move onto the next exon for further matching
            nextRegionStart = transData.exons().get(i + 1).Start;
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
        if(readRegions.isEmpty())
            return false;

        int regionsStart = readRegions.get(0)[SE_START];
        int regionsEnd = readRegions.get(readRegions.size() - 1)[SE_END];

        if(!positionsOverlap(regionsStart, regionsEnd, transData.TransStart, transData.TransEnd))
            return false;

        if(requiredMatchType == SHORT)
        {
            // region must lie within an exon
            return transData.exons().stream().anyMatch(x -> positionsWithin(regionsStart, regionsEnd, x.Start, x.End));
        }
        else
        {
            // first check for matching splice junctions
            int exonIndex = 0;
            for(int[] spliceJunction : spliceJunctions)
            {
                int spliceStart = spliceJunction[SE_START];
                int spliceEnd = spliceJunction[SE_END];
                boolean matched = false;

                for (; exonIndex < transData.exons().size() - 1; ++exonIndex)
                {
                    ExonData exon = transData.exons().get(exonIndex);
                    ExonData nextExon = transData.exons().get(exonIndex + 1);

                    if (exon.End == spliceStart && nextExon.Start == spliceEnd)
                    {
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    return false;
            }

            // now check the each of the read regions falls within an exon - exons can be skipped
            // and whether the total implief fragment length from the POV of this transcript makes it too long
            exonIndex = 0;
            for(int[] readRegion : readRegions)
            {
                int regionStart = readRegion[SE_START];
                int regionEnd = readRegion[SE_END];
                boolean matched = false;

                for(; exonIndex < transData.exons().size(); ++exonIndex)
                {
                    final ExonData exon = transData.exons().get(exonIndex);

                    if(regionStart > exon.End) // keep searching
                        continue;

                    if(regionEnd < exon.Start) // would be intronic for this transcript
                        return false;

                    if(positionsOverlap(regionStart, regionEnd, exon.Start, exon.End))
                    {
                        if(!positionsWithin(regionStart, regionEnd, exon.Start, exon.End))
                            return false; // breaches the exon boundaries

                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    return false;
            }

            int impliedFragLength = calcFragmentLength(transData, regionsStart, regionsEnd);
            if(impliedFragLength > mConfig.MaxFragmentLength)
                return false;
        }

        return true;
    }

    public static void formTranscriptDefinitions(final List<CategoryCountsData> categoryCountsData, ExpectedRatesData expRatesData)
    {
        // convert fragment counts in each category per transcript into the equivalent of a signature per transcript
        collectCategories(categoryCountsData, expRatesData);

        final Map<String,List<CategoryCountsData>> transGeneCountsMap = createTransComboDataMap(categoryCountsData);

        int categoryCount = expRatesData.Categories.size();

        transGeneCountsMap.keySet().forEach(x -> expRatesData.TranscriptIds.add(x));

        expRatesData.initialiseTranscriptDefinitions();

        for(int transIndex = 0; transIndex < expRatesData.TranscriptIds.size(); ++transIndex)
        {
            final String transId = expRatesData.TranscriptIds.get(transIndex);

            double[] categoryCounts = new double[categoryCount];

            final List<CategoryCountsData> transCounts = transGeneCountsMap.get(transId);

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

    private static void collectCategories(
            final List<CategoryCountsData> categoryCountsData, ExpectedRatesData expRatesData)
    {
        for(final CategoryCountsData tcData : categoryCountsData)
        {
            final String transKey = tcData.combinedKey();

            if(tcData.fragmentCount() > 0 || tcData.transcriptIds().isEmpty()) // force inclusion of unspliced gene categories
            {
                expRatesData.addCategory(transKey);
            }
        }
    }

    public static final Map<String,List<CategoryCountsData>> createTransComboDataMap(final List<CategoryCountsData> categoryCountsData)
    {
        final Map<String,List<CategoryCountsData>> transGeneCountsMap = Maps.newHashMap();

        for(final CategoryCountsData tcData : categoryCountsData)
        {
            final Set<String> geneTransNames = Sets.newHashSet();
            tcData.unsplicedGeneIds().forEach(x -> geneTransNames.add(x));
            tcData.transcriptIds().forEach(x -> geneTransNames.add(String.valueOf(x)));

            for(String geneTransName : geneTransNames)
            {
                List<CategoryCountsData> countsList = transGeneCountsMap.get(geneTransName);
                if(countsList == null)
                {
                    transGeneCountsMap.put(geneTransName, Lists.newArrayList(tcData));
                }
                else
                {
                    countsList.add(tcData);
                }
            }
        }

        return transGeneCountsMap;
    }

    public void setFragmentLengthData(int length, int frequency)
    {
        mCurrentFragSize = length;
        mCurrentFragFrequency = frequency;
        mReadLength = mConfig.ReadLength;
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
                writer.write("GeneSetId,Category");

                for(FragmentSize fragLength : config.FragmentSizeData)
                {
                    writer.write(String.format(",%s%d", EXP_COUNT_LENGTH_HEADER, fragLength.Length));
                }
            }
            else
            {
                writer.write("GeneSetId,Category,Rate");
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

    private synchronized static void writeExpectedCounts(
            final BufferedWriter writer, final String collectionId, final List<CategoryCountsData> categoryCounts)
    {
        if(writer == null)
            return;

        try
        {
            for(CategoryCountsData tcData : categoryCounts)
            {
                final double[] lengthCounts = tcData.fragmentCountsByLength();

                writer.write(String.format("%s,%s,%.0f", collectionId, tcData.combinedKey(), lengthCounts[0]));

                for (int i = 1; i < lengthCounts.length; ++i)
                {
                    writer.write(String.format(",%.0f", lengthCounts[i]));
                }

                writer.newLine();
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
            final Matrix transcriptDefinitions = expExpData.getTranscriptDefinitions();
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

    @VisibleForTesting
    public Map<String,List<CategoryCountsData>> getTransComboDataMap() { return mTransCategoryCountsMap; }

}
