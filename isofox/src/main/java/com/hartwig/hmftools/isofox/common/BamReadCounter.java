package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.ChromosomeTaskExecutor.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

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
    private final ResultsWriter mResultsWriter;

    private final int[] mCurrentGenesRange;
    private int mTotalReadCount;
    private int mCurrentGeneReadCount;
    private FragmentTypeCounts mFragmentTypeCounts;
    private int mSecondaryReads;
    private String mChromosome;
    private final List<GeneData> mGeneDataList;
    private String mCurrentGenes;
    private final int[] mMaqQualFrequencies;

    public BamReadCounter(final IsofoxConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mResultsWriter = resultsWriter;

        mGeneDataList = Lists.newArrayList();
        mChromosome = "";

        mCurrentGenesRange = new int[SE_PAIR];
        mTotalReadCount = 0;
        mCurrentGeneReadCount = 0;
        mCurrentGenes = "";
        mSecondaryReads = 0;
        mFragmentTypeCounts = new FragmentTypeCounts();
        mMaqQualFrequencies = new int[4];
    }

    public void initialise(final String chromosome, final List<GeneData> geneDataList)
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
        if(mGeneDataList.isEmpty())
            return;

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.info("processing reads for chromosome({}) geneCount({})", mChromosome, mGeneDataList.size());

        final List<GeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < mGeneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(mGeneDataList, currentGeneIndex, overlappingGenes);

            // if(overlappingGenes.stream().anyMatch(x -> mConfig.Filters.EnrichedGeneIds.contains(x.GeneId)))
            //    continue;

            mCurrentGenesRange[SE_START] = 0;
            mCurrentGenesRange[SE_END] = 0;

            for(int i = 0; i < overlappingGenes.size(); ++i)
            {
                GeneData geneData = overlappingGenes.get(i);

                mCurrentGenesRange[SE_START] = i == 0 ? geneData.GeneStart : min(geneData.GeneStart, mCurrentGenesRange[SE_START]);
                mCurrentGenesRange[SE_END] = i == 0 ? geneData.GeneEnd : max(geneData.GeneEnd, mCurrentGenesRange[SE_END]);
            }

            if(currentGeneIndex >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chromosome({}) processed {} genes, totalReads({})",
                        mChromosome, currentGeneIndex, mTotalReadCount);
            }

            mCurrentGenes = overlappingGenes.get(0).GeneId;
            mCurrentGeneReadCount = 0;

            final List<ChrBaseRegion> regions = Lists.newArrayList(new ChrBaseRegion(mChromosome, mCurrentGenesRange));

            mBamSlicer.slice(mSamReader, regions, this::processBamRead);
        }

        ISF_LOGGER.info("chromosome({}) processing complete: total({}) duplicates({}) chimeric({}) secondaries({}) mapQuals(0={} 1={} 2={} 3={})",
                mChromosome, mTotalReadCount, mFragmentTypeCounts.typeCount(DUPLICATE), mFragmentTypeCounts.typeCount(CHIMERIC), mSecondaryReads,
                mMaqQualFrequencies[0], mMaqQualFrequencies[1], mMaqQualFrequencies[2], mMaqQualFrequencies[3]);
    }

    private void processBamRead(final SAMRecord record)
    {
        if(!positionWithin(record.getStart(), mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]))
            return;

        if(GeneRegionFilters.inExcludedRegion(mConfig.Filters.ExcludedRegion, record))
            return;

        ++mTotalReadCount;
        ++mCurrentGeneReadCount;
        mFragmentTypeCounts.addCount(TOTAL);

        if(record.getDuplicateReadFlag())
            mFragmentTypeCounts.addCount(DUPLICATE);

        if((record.getFlags() & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0)
            mFragmentTypeCounts.addCount(CHIMERIC);

        if(record.isSecondaryAlignment())
            ++mSecondaryReads;

        if(record.getMappingQuality() <= 3)
        {
            mMaqQualFrequencies[record.getMappingQuality()]++;
        }

        if(mConfig.GeneReadLimit > 0 && mCurrentGeneReadCount > mConfig.GeneReadLimit)
        {
            mBamSlicer.haltProcessing();
            ISF_LOGGER.info("chromosome({}) gene({}) halting processing after {} reads", mChromosome, mCurrentGenes, mCurrentGeneReadCount);
        }

        if(mConfig.WriteReadData)
            writeReadData(mResultsWriter.getReadDataWriter(), record, mCurrentGenes);
    }

    public static BufferedWriter createReadDataWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("read_data.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,ReadId,Chromosome,PosStart,PosEnd,Cigar,Flags,InsertSize");
            writer.write(",MateChr,MatePosStart,FirstInPair,ReadReversed,Duplicate,Secondary,Supplementary,SuppData");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeReadData( final BufferedWriter writer, final SAMRecord record, final String geneId)
    {
        try
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            writer.write(String.format("%s,%s,%s,%d,%d,%s,%d,%d,%s,%d",
                    geneId, record.getReadName(), record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd(),
                    record.getCigarString(), record.getFlags(), record.getInferredInsertSize(),
                    record.getMateReferenceName(), record.getMateAlignmentStart()));

            writer.write(String.format(",%s,%s,%s,%s,%s,%s",
                    record.getFirstOfPairFlag(), record.getReadNegativeStrandFlag(), record.getDuplicateReadFlag(),
                    record.getSecondOfPairFlag(), record.getSupplementaryAlignmentFlag(), suppData != null ? suppData.asCsv() : "N/A"));

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write read data file: {}", e.toString());
        }
    }
}
