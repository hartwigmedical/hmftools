package com.hartwig.hmftools.isofox.refdata;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.FragmentAllocator.calcFragmentLength;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.LONG;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.UNSPLICED;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.formTranscriptDefinitions;
import static com.hartwig.hmftools.isofox.refdata.RefDataWriter.writeExpectedCounts;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.expression.CategoryCountsData;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesData;

public class ExpectedCountsGenerator
{
    private final RefDataConfig mConfig;

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

    public ExpectedCountsGenerator(final RefDataConfig config, final RefDataWriter resultsWriter)
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

            for(TranscriptData transData : transDataList)
            {
                boolean endOfTrans = false;
                List<TranscriptData> candidateTrans = Lists.newArrayList(transDataList);
                candidateTrans.remove(transData);

                for(ExonData exon : transData.exons())
                {
                    for(int startPos = exon.Start; startPos <= exon.End; ++startPos)
                    {
                        cullTranscripts(candidateTrans, startPos);

                        if(!allocateTranscriptCounts(transData, transDataList, startPos))
                        {
                            endOfTrans = true;
                            break;
                        }
                    }

                    if(endOfTrans)
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

                for(int startPos = regionStart; startPos <= regionEnd - mCurrentFragSize; ++startPos)
                {
                    final List<String> unsplicedGenes = findUnsplicedGenes(startPos);

                    // cull the set of possible transcripts
                    cullTranscripts(candidateTrans, startPos);

                    if(startPos <= currentExonicEnd)
                    {
                        // check possible transcript exonic matches
                        allocateUnsplicedCounts(transDataList, startPos, unsplicedGenes);
                    }
                    else
                    {
                        // check for purely intronic fragments
                        if(startPos < nextExonicStart)
                        {
                            addUnsplicedCountsData(emptyTrans, unsplicedGenes);
                        }
                        else
                        {
                            ++exonicRegionIndex;
                            currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];

                            if(exonicRegionIndex < commonExonicRegions.size() - 1)
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

        writeExpectedCounts(mExpRateWriter, geneCollection.chrId(), mTransCategoryCounts);
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

            matchingCounts.initialiseLengthCounts(mConfig.FragmentSizeData.size());

            transComboDataList.add(matchingCounts);
        }

        matchingCounts.addFragLengthCounts(mCurrentFragFrequency, mFragSizeIndex);
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

                if(readsAdded == 2)
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

                for(; exonIndex < transData.exons().size() - 1; ++exonIndex)
                {
                    ExonData exon = transData.exons().get(exonIndex);
                    ExonData nextExon = transData.exons().get(exonIndex + 1);

                    if(exon.End == spliceStart && nextExon.Start == spliceEnd)
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

    @VisibleForTesting
    public void setFragmentLengthData(int length, int frequency)
    {
        mCurrentFragSize = length;
        mCurrentFragFrequency = frequency;
        mReadLength = mConfig.ReadLength;
    }

    @VisibleForTesting
    public void populateExpectedRates()
    {
        formTranscriptDefinitions(mTransCategoryCounts, mCurrentExpRatesData);
    }

    @VisibleForTesting
    public Map<String,List<CategoryCountsData>> getTransComboDataMap() { return mTransCategoryCountsMap; }
}
