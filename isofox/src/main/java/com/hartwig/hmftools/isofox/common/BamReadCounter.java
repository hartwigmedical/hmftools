package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import static com.hartwig.hmftools.isofox.BamFragmentReader.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// simple BAM read counter, used for experimental purposes at the moment
public class BamReadCounter implements Callable
{
    private final IsofoxConfig mConfig;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final int[] mCurrentGenesRange;
    private int mTotalReadCount;
    private int mCurrentGeneReadCount;
    private int[] mReadTypeCounts;
    private int mSecondaryReads;
    private String mChromosome;
    private final List<EnsemblGeneData> mGeneDataList;
    private String mCurrentGenes;
    private final int[] mMaqQualFrequencies;

    public BamReadCounter(final IsofoxConfig config)
    {
        mConfig = config;
        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);

        mGeneDataList = Lists.newArrayList();
        mChromosome = "";

        mCurrentGenesRange = new int[SE_PAIR];
        mTotalReadCount = 0;
        mCurrentGeneReadCount = 0;
        mCurrentGenes = "";
        mSecondaryReads = 0;
        mReadTypeCounts = new int[typeAsInt(FragmentType.MAX)];
        mMaqQualFrequencies = new int[4];
    }

    public void initialise(final String chromosome, final List<EnsemblGeneData> geneDataList)
    {
        mChromosome = chromosome;
        mGeneDataList.clear();
        mGeneDataList.addAll(geneDataList);
    }

    @Override
    public Long call()
    {
        processBam();
        return (long)0;
    }

    private void processBam()
    {
        if (mGeneDataList.isEmpty())
            return;

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.info("processing reads for chromosome({}) geneCount({})", mChromosome, mGeneDataList.size());

        final List<EnsemblGeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < mGeneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(mGeneDataList, currentGeneIndex, overlappingGenes);

            if(overlappingGenes.stream().anyMatch(x -> mConfig.EnrichedGeneIds.contains(x.GeneId)))
                continue;

            mCurrentGenesRange[SE_START] = 0;
            mCurrentGenesRange[SE_END] = 0;

            for (int i = 0; i < overlappingGenes.size(); ++i)
            {
                EnsemblGeneData geneData = overlappingGenes.get(i);

                mCurrentGenesRange[SE_START] = i == 0 ? geneData.GeneStart : min(geneData.GeneStart, mCurrentGenesRange[SE_START]);
                mCurrentGenesRange[SE_END] = i == 0 ? geneData.GeneEnd : max(geneData.GeneEnd, mCurrentGenesRange[SE_END]);
            }

            if(currentGeneIndex >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.info("chromosome({}) processed {} genes, totalReads({})",
                        mChromosome, currentGeneIndex, mTotalReadCount);
            }

            mCurrentGenes = overlappingGenes.get(0).GeneId;
            mCurrentGeneReadCount = 0;

            final List<BaseRegion> regions = Lists.newArrayList(new BaseRegion(mChromosome, mCurrentGenesRange));

            mBamSlicer.slice(mSamReader, regions, this::processBamRead);
        }

        ISF_LOGGER.info("chromosome({}) processing complete: total({}) duplicates({}) chimeric({}) secondaries({}) mapQuals(0={} 1={} 2={} 3={})",
                mChromosome, mTotalReadCount, mReadTypeCounts[typeAsInt(DUPLICATE)], mReadTypeCounts[typeAsInt(CHIMERIC)], mSecondaryReads,
                mMaqQualFrequencies[0], mMaqQualFrequencies[1], mMaqQualFrequencies[2], mMaqQualFrequencies[3]);
    }

    private void processBamRead(@NotNull final SAMRecord read)
    {
        ++mTotalReadCount;
        ++mCurrentGeneReadCount;
        ++mReadTypeCounts[typeAsInt(TOTAL)];

        if(read.getDuplicateReadFlag())
            ++mReadTypeCounts[typeAsInt(DUPLICATE)];

        if((read.getFlags() & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0)
            ++mReadTypeCounts[typeAsInt(CHIMERIC)];

        if(read.isSecondaryAlignment())
            ++mSecondaryReads;

        if(read.getMappingQuality() <= 3)
        {
            mMaqQualFrequencies[read.getMappingQuality()]++;
        }

        if(mConfig.GeneReadLimit > 0 && mCurrentGeneReadCount > mConfig.GeneReadLimit)
        {
            mBamSlicer.haltProcessing();
            ISF_LOGGER.info("chromosome({}) gene({}) halting processing after {} reads", mChromosome, mCurrentGenes, mCurrentGeneReadCount);
        }
    }

}
