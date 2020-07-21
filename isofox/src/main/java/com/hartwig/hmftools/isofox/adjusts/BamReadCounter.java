package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.ChromosomeGeneTask.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.FragmentType;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// simple BAM read counter, used for experimental purposes at the moment

public class BamReadCounter
{
    private final IsofoxConfig mConfig;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final int[] mCurrentGenesRange;
    private int mTotalReadCount;
    private int mCurrentGeneReadCount;
    private int[] mReadTypeCounts;
    private String mChromosome;
    private String mCurrentGenes;
    private final FragmentTracker mFragmentTracker;

    public BamReadCounter(final IsofoxConfig config)
    {
        mConfig = config;
        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false, false);

        mCurrentGenesRange = new int[SE_PAIR];
        mTotalReadCount = 0;
        mCurrentGeneReadCount = 0;
        mCurrentGenes = "";
        mChromosome = "";
        mReadTypeCounts = new int[typeAsInt(FragmentType.MAX)];
        mFragmentTracker = new FragmentTracker();
    }

    public void processBam(final String chromosome, final List<EnsemblGeneData> geneDataList)
    {
        if (geneDataList == null || geneDataList.isEmpty())
            return;

        mChromosome = chromosome;

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.info("processing reads for chromosome({}) geneCount({})", chromosome, geneDataList.size());

        final List<EnsemblGeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < geneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(geneDataList, currentGeneIndex, overlappingGenes);

            if(overlappingGenes.stream().anyMatch(x -> mConfig.EnrichedGeneIds.contains(x.GeneId)))
                continue;

            mFragmentTracker.clear();
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
                        chromosome, currentGeneIndex, mTotalReadCount);
            }

            mCurrentGenes = overlappingGenes.get(0).GeneId;
            mCurrentGeneReadCount = 0;

            List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]));

            mBamSlicer.slice(mSamReader, regions, this::processBamRead);
        }

        ISF_LOGGER.info("chromosome({}) processing complete: total({}) duplicates({}) chimeric({})",
                chromosome, mTotalReadCount, mReadTypeCounts[typeAsInt(DUPLICATE)], mReadTypeCounts[typeAsInt(CHIMERIC)]);
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

        if(mConfig.GeneReadLimit > 0 && mCurrentGeneReadCount > mConfig.GeneReadLimit)
        {
            mBamSlicer.haltProcessing();
            ISF_LOGGER.info("chromosome({}) gene({}) halting processing after {} reads", mChromosome, mCurrentGenes, mCurrentGeneReadCount);
        }
    }

}
