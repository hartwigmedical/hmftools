package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_MAX;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_UNSPLICED;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExpectedExpressionRates
{
    private final RnaExpConfig mConfig;

    private final int[] mIntronicCounts;
    private final Map<String,List<TranscriptComboData>> mTransComboData;

    private static final Logger LOGGER = LogManager.getLogger(ExpectedExpressionRates.class);

    public ExpectedExpressionRates(final RnaExpConfig config)
    {
        mConfig = config;

        mTransComboData = Maps.newHashMap();
        mIntronicCounts = new int[TC_MAX];
    }

    public void generate(final GeneReadData geneReadData)
    {
        long regionStart = geneReadData.getTranscriptsRange()[SE_START];
        long regionEnd = geneReadData.getTranscriptsRange()[SE_END];

        int fragmentLength = mConfig.MedianFragmentLength;

        final List<long[]> commonExonicRegions = geneReadData.getCommonExonicRegions();

        int exonicRegionIndex = 0;
        long currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];

        for(long startPos = regionStart; startPos < regionEnd - fragmentLength; ++startPos)
        {
            if(startPos <= currentExonicEnd)
            {
                // check transcript allocations
                allocateTranscriptCounts(geneReadData, startPos);
                continue;
            }

            if(exonicRegionIndex < commonExonicRegions.size() - 1)
            {
                long nextExonicStart = commonExonicRegions.get(exonicRegionIndex + 1)[SE_START];

                if (startPos < nextExonicStart)
                {
                    ++mIntronicCounts[TC_UNSPLICED];
                }
                else
                {
                    ++exonicRegionIndex;
                    currentExonicEnd = commonExonicRegions.get(exonicRegionIndex)[SE_END];

                    // check transcript allocations
                    allocateTranscriptCounts(geneReadData, startPos);
                }
            }
        }
    }

    private void allocateTranscriptCounts(final GeneReadData geneReadData, long startPos)
    {
        final List<TranscriptData> transDataList = geneReadData.getTranscripts();

        boolean matchedTrans = false;

        for(TranscriptData transData : transDataList)
        {
            List<long[]> readRegions = Lists.newArrayList();

            boolean spanExons = generateImpliedFragment(transData, startPos, readRegions);

            if(readRegions.isEmpty())
                continue;

            matchedTrans = true;

            final List<String> longAndSplicedTrans = Lists.newArrayList();
            final List<String> shortTrans = Lists.newArrayList();

            if(spanExons)
                longAndSplicedTrans.add(transData.TransName);
            else
                shortTrans.add(transData.TransName);

            // now check whether these regions are supported by each other transcript
            for(TranscriptData otherTransData : transDataList)
            {
                if(transData == otherTransData)
                    continue;

                if(readsSupportFragment(otherTransData, readRegions, spanExons))
                {
                    if(spanExons)
                        longAndSplicedTrans.add(transData.TransName);
                    else
                        shortTrans.add(transData.TransName);
                }
            }

            if(!longAndSplicedTrans.isEmpty())
            {
                addTransComboData(transData.TransName, longAndSplicedTrans, TC_SPLICED);
            }
            else
            {
                addTransComboData(transData.TransName, shortTrans, TC_SHORT);
            }
        }

        if(!matchedTrans)
        {
            ++mIntronicCounts[TC_UNSPLICED];
        }
    }

    private void addTransComboData(final String transId, final List<String> transcripts, int transMatchType)
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

        ++matchingCounts.getCounts()[transMatchType];
    }

    public boolean generateImpliedFragment(final TranscriptData transData, long startPos, List<long[]> readRegions)
    {
        readRegions.clear();

        // set out the fragment reads either within a single exon or spanning one or more
        int fragmentLength = mConfig.MedianFragmentLength;

        int exonCount = transData.exons().size();
        final ExonData lastExon = transData.exons().get(exonCount - 1);

        if(startPos + fragmentLength - 1 > lastExon.ExonEnd)
            return false;

        int readLength = mConfig.ReadLength;

        int remainingReadBases = readLength;
        int remainingInterimBases = fragmentLength - 2 * readLength + 1;
        long nextRegionStart = startPos;
        int readsAdded = 0;
        boolean spansExons = false;

        for(int i = 0; i < exonCount; ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(nextRegionStart > exon.ExonEnd)
                continue;

            if(!readRegions.isEmpty())
                spansExons = true;

            if(readsAdded == 1 && remainingInterimBases > 0)
            {
                if(nextRegionStart + remainingInterimBases - 1 >= exon.ExonEnd)
                {
                    if(i >= exonCount - 1)
                    {
                        readRegions.clear();
                        return false;
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

            if(remainingReadBases == 0)
            {
                ++readsAdded;

                if (readsAdded == 2)
                    break;

                remainingReadBases = readLength;
            }

            // is the remainer of this exon long enough to match again?
            nextRegionStart = regionEnd + remainingInterimBases;

            if(regionEnd == exon.ExonEnd || nextRegionStart > exon.ExonEnd)
            {
                if(i == exonCount - 1)
                {
                    readRegions.clear();
                    return false;
                }

                // will move onto the next exon for further matching
                nextRegionStart = transData.exons().get(i + 1).ExonStart;

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
            readRegions.add(new long[] {nextRegionStart, regionEnd});

            if(remainingReadBases == 0)
                break;

            if(i == exonCount - 1)
            {
                readRegions.clear();
                return false;
            }

            // will move onto the next exon for further matching
            nextRegionStart = transData.exons().get(i + 1).ExonStart;
        }

        return spansExons;
    }

    public static boolean readsSupportFragment(final TranscriptData transData, List<long[]> readRegions, boolean requireSpan)
    {
        long regionsStart = readRegions.get(0)[SE_START];
        long regionsEnd = readRegions.get(readRegions.size() - 1)[SE_END];

        if(!requireSpan)
        {
            return (transData.exons().stream().anyMatch(x -> x.ExonStart <= regionsStart && x.ExonEnd >= regionsEnd));
        }
        else
        {
            // none of the reads can breach an exon boundary
            int straddlingExons = 0;
            for(int i = 0; i < readRegions.size(); ++i)
            {
                long regionStart = readRegions.get(i)[SE_START];
                long regionEnd = readRegions.get(i)[SE_END];

                boolean withinExon = transData.exons().stream().anyMatch(x -> x.ExonStart <= regionStart && x.ExonEnd >= regionEnd);

                if(!withinExon)
                    return false;
            }

            // cannot be the same exon
            if(transData.exons().stream().anyMatch(x -> x.ExonStart <= regionsStart && x.ExonEnd >= regionsEnd))
                return false;
        }

        return true;
    }

}
