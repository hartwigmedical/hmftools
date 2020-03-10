package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// calculate fragment length distribution for a sample
// this can be done either independently from fragment counting or lengths can be registering during that process

public class FragmentSizeCalcs
{
    private final RnaExpConfig mConfig;
    private final SvGeneTranscriptCollection mGeneTransCache;

    private final List<int[]> mFragmentLengths;
    private final List<int[]> mFragmentLengthsByGene;
    private int mMaxReadLength;

    private BufferedWriter mFragLengthWriter;
    private BufferedWriter mFragReadWriter;

    private EnsemblGeneData mCurrentGeneData;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentReadCount;
    private int mTotalReadCount;
    private int mProcessedFragments;

    private PerformanceCounter mPerfCounter;

    public FragmentSizeCalcs(final RnaExpConfig config, final SvGeneTranscriptCollection geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mCurrentGeneData = null;
        mCurrentTransDataList = null;
        mCurrentReadCount = 0;
        mTotalReadCount = 0;
        mProcessedFragments = 0;
        mMaxReadLength = 0;

        mFragLengthWriter = writer;
        mFragReadWriter = null;

        mFragmentLengths = Lists.newArrayList();
        mFragmentLengthsByGene = Lists.newArrayList();
        mPerfCounter = new PerformanceCounter("FragLengthDist");
    }

    public void close()
    {
        closeBufferedWriter(mFragReadWriter);

        if(RE_LOGGER.isDebugEnabled())
            mPerfCounter.logStats();
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }
    public final int getMaxReadLength() { return mMaxReadLength; }

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 20;
    private static final int MAX_TRAN_EXONS = 20;
    private static final int MAX_GENE_READ_COUNT = 1000; // to avoid impact of highly enriched genes

    private static int FRAG_LENGTH_CAP = 3000; // to prevent map blowing out in size

    public void calcSampleFragmentSize(final String chromosome, final List<EnsemblGeneData> geneDataList)
    {
        if (geneDataList == null || geneDataList.isEmpty())
            return;

        // measure fragment lengths for non-split reads in purely intronic regions
        if(mGeneTransCache == null)
        {
            RE_LOGGER.error("fragment length calculator missing gene cache");
            return;
        }

        if (mConfig.FragmentLengthMinCount > 0 && mProcessedFragments >= mConfig.FragmentLengthMinCount)
            return;

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        BamSlicer bamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        // walk through each chromosome, ignoring any gene which overlaps the previous gene
        RE_LOGGER.debug("calculating fragment size for chromosome({}) geneCount({})", chromosome, geneDataList.size());

        long lastGeneEnd = 0;

        for (int i = 0; i < geneDataList.size(); ++i)
        {
            EnsemblGeneData geneData = geneDataList.get(i);

            if(mConfig.ExcludedGeneIds.contains(geneData.GeneId))
                continue;

            if (geneData.GeneStart < lastGeneEnd)
                continue;

            long geneLength = geneData.GeneEnd - geneData.GeneStart;

            if (geneLength < MIN_GENE_LENGTH || geneLength > MAX_GENE_LENGTH)
                continue;

            mCurrentTransDataList = mGeneTransCache.getTranscripts(geneData.GeneId).stream()
                    .filter(x -> x.exons().size() <= MAX_TRAN_EXONS)
                    .collect(Collectors.toList());

            if (mCurrentTransDataList.isEmpty() || mCurrentTransDataList.size() > MAX_GENE_TRANS)
                continue;

            if(i > 0 && (i % 100) == 0)
            {
                RE_LOGGER.debug("chromosome({}) processed {} genes, lastGenePos({}) fragCount({}) totalReads({})",
                        chromosome, i, lastGeneEnd, mProcessedFragments, mTotalReadCount);
            }

            lastGeneEnd = geneData.GeneEnd;

            mPerfCounter.start();

            mCurrentReadCount = 0;
            mCurrentGeneData = geneData;
            List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, geneData.GeneStart, geneData.GeneEnd));

            bamSlicer.slice(samReader, regions, this::processBamRead);

            mPerfCounter.stop();

            if(mConfig.FragmentLengthsByGene)
            {
                writeFragmentLengths(mFragLengthWriter, mFragmentLengthsByGene, geneData);
                mFragmentLengthsByGene.clear();
            }

            if (mConfig.FragmentLengthMinCount > 0 && mProcessedFragments >= mConfig.FragmentLengthMinCount)
            {
                RE_LOGGER.debug("max fragment length samples reached: {}", mProcessedFragments);
                break;
            }
        }
    }

    private void processBamRead(@NotNull final SAMRecord record)
    {
        ++mCurrentReadCount;
        ++mTotalReadCount;

        if(mTotalReadCount > 0 && (mTotalReadCount % 10000) == 0)
        {
            RE_LOGGER.trace("currentGene({}:{}) totalReads({})", mCurrentGeneData.GeneId, mCurrentGeneData.GeneName, mTotalReadCount);
        }

        if(mCurrentReadCount >= MAX_GENE_READ_COUNT)
            return;

        if(!isCandidateRecord(record))
            return;

        mMaxReadLength = max(mMaxReadLength, record.getReadLength());

        long posStart = record.getStart();
        long posEnd = record.getEnd();

        for(final TranscriptData transData : mCurrentTransDataList)
        {
            if(transData.exons().stream().anyMatch(x -> !(posStart > x.ExonEnd || posEnd < x.ExonStart)))
                return;
        }

        addFragmentLength(record, mFragmentLengths);

        if(mConfig.FragmentLengthsByGene)
            addFragmentLength(record, mFragmentLengthsByGene);

        if(mConfig.WriteFragmentReads)
            writeReadData(mCurrentGeneData, record);
    }

    private boolean isCandidateRecord(final SAMRecord record)
    {
        if(!record.getFirstOfPairFlag())
            return false;

        // ignore translocations and inversions
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
            return false;

        // ignore split or unmapped reads
        if(record.getCigar() == null || record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return false;

        return true;
    }

    private synchronized void addFragmentLength(final SAMRecord record, final List<int[]> fragmentLengths)
    {
        int fragmentLength = min(abs(record.getInferredInsertSize()), FRAG_LENGTH_CAP);

        if(fragmentLength == 0)
            return;

        int index = 0;
        boolean exists = false;
        while(index < fragmentLengths.size())
        {
            final int[] fragLengthCount = fragmentLengths.get(index);

            if(fragLengthCount[FL_LENGTH] < fragmentLength)
            {
                ++index;
                continue;
            }

            if(fragLengthCount[FL_LENGTH] == fragmentLength)
            {
                ++fragLengthCount[FL_FREQUENCY];
                exists = true;
            }

            break;
        }

        if(!exists)
        {
            int[] newFragLengthCount = { fragmentLength, 1 };
            fragmentLengths.add(index, newFragLengthCount);
        }

        ++mProcessedFragments;
    }

    public void setConfigFragmentLengthData()
    {
        if(mMaxReadLength > 0)
        {
            mConfig.ReadLength = mMaxReadLength;
            RE_LOGGER.info("max read length({}) set", mMaxReadLength);
        }
        else
        {
            RE_LOGGER.warn("max read length not determined from fragment length calcs");
        }

        final List<int[]> lengthFrequencies = mConfig.ExpRateFragmentLengths;

        int currentRangeMin = 0;
        int currentRangeMax = 0;


        for(int i = 0; i < lengthFrequencies.size(); ++i)
        {
            int[] lengthFrequency = lengthFrequencies.get(i);

            currentRangeMin = (i == 0) ? 0 : currentRangeMax + 1;

            if(i == lengthFrequencies.size() - 1)
            {
                currentRangeMax = mConfig.MaxFragmentLength - 1;
            }
            else
            {
                int[] nextLengthFrequency = lengthFrequencies.get(i + 1);
                currentRangeMax = (lengthFrequency[FL_LENGTH] + nextLengthFrequency[FL_LENGTH]) / 2;
            }

            int lengthCount = 0;

            for (final int[] fragLengthCount : mFragmentLengths)
            {
                if(fragLengthCount[FL_LENGTH] >= currentRangeMin && fragLengthCount[FL_LENGTH] <= currentRangeMax)
                {
                    lengthCount += fragLengthCount[FL_FREQUENCY];
                }
            }

            lengthFrequency[FL_FREQUENCY] = lengthCount;
            RE_LOGGER.info("fragmentLength({}) frequency({})", lengthFrequency[FL_LENGTH], lengthCount);
        }
    }

    public static BufferedWriter createFragmentLengthWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.FragmentLengthsByGene ?
                    config.formOutputFile("frag_length_by_gene.csv") : config.formOutputFile("frag_length.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if (config.FragmentLengthsByGene)
            {
                writer.write("GeneId,GeneName,Chromosome,");
            }

            writer.write("FragmentLength,Count");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write fragment length file: {}", e.toString());
            return null;
        }
    }

    public void writeFragmentLengths(final EnsemblGeneData geneData)
    {
        writeFragmentLengths(mFragLengthWriter, mFragmentLengths, null);
    }

    private synchronized static void writeFragmentLengths(
            final BufferedWriter writer, @Nullable final List<int[]> fragmentLengths, final EnsemblGeneData geneData)
    {
        if(writer == null)
            return;

        if(fragmentLengths.isEmpty())
            return;

        try
        {
            for (final int[] fragLengthCount : fragmentLengths)
            {
                if(geneData != null)
                {
                    writer.write(String.format("%s,%s,%s,",
                            geneData.GeneId, geneData.GeneName, geneData.Chromosome));
                }

                writer.write(String.format("%d,%d", fragLengthCount[FL_LENGTH], fragLengthCount[FL_FREQUENCY]));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }

    }

    private void writeReadData(final EnsemblGeneData geneData, final SAMRecord read)
    {
        if(mConfig.Threads > 1) // not currently synchronised
            return;

        try
        {
            if(mFragReadWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("frag_length_reads.csv");

                mFragReadWriter = createBufferedWriter(outputFileName, false);
                mFragReadWriter.write("GeneId,GeneName,Chromosome,ReadId,PosStart,PosEnd,Cigar");
                mFragReadWriter.write(",InsertSize,ReadLength,OtherReadStart,HasReadOverlap,ReadBases");
                mFragReadWriter.newLine();
            }

            mFragReadWriter.write(String.format("%s,%s,%s",
                    geneData.GeneId, geneData.GeneName, geneData.Chromosome));

            long otherReadStart = read.getMateAlignmentStart();
            int readLength = read.getReadLength();

            boolean hasReadOverlap;
            if(otherReadStart < read.getStart())
            {
                hasReadOverlap = otherReadStart + readLength - 1 >= read.getStart();
            }
            else
            {
                hasReadOverlap = read.getEnd() >= otherReadStart;
            }

            mFragReadWriter.write(String.format(",%s,%d,%d,%s",
                    read.getReadName(), read.getStart(), read.getEnd(), read.getCigar().toString()));

            mFragReadWriter.write(String.format(",%d,%d,%d,%s,%s",
                    read.getInferredInsertSize(), readLength, otherReadStart, hasReadOverlap, read.getReadString()));

            mFragReadWriter.newLine();
        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write fragment length read data file: {}", e.toString());
        }
    }

}
