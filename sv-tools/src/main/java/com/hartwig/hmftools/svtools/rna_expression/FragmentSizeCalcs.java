package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class FragmentSizeCalcs
{
    private final RnaExpConfig mConfig;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final RnaBamReader mRnaBamReader;

    private final List<int[]> mFragmentLengths;
    private BufferedWriter mWriter;

    // indices into fragment length frequency array
    public static final int FL_LENGTH = 0;
    public static final int FL_FREQUENCY = 1;

    private EnsemblGeneData mCurrentGeneData;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentReadCount;
    private int mTotalReadCount;
    private int mProcessedFragments;

    private static final Logger LOGGER = LogManager.getLogger(FragmentSizeCalcs.class);

    public FragmentSizeCalcs(final RnaExpConfig config, final SvGeneTranscriptCollection geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;
        mRnaBamReader = new RnaBamReader(config);

        mCurrentGeneData = null;
        mCurrentTransDataList = null;
        mCurrentReadCount = 0;
        mTotalReadCount = 0;
        mProcessedFragments = 0;
        mWriter = null;

        mFragmentLengths = Lists.newArrayList();
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }
    public boolean enabled() { return mConfig.WriteFragmentLengths || mConfig.FragmentLengthMinCount > 0; }

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 20;
    private static final int MAX_TRAN_EXONS = 20;
    private static final int MAX_GENE_READ_COUNT = 1000; // to avoid impact of highly enriched genes

    private static int FRAG_LENGTH_CAP = 5000; // to prevent map blowing out in size

    public void calcSampleFragmentSize()
    {
        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> geneDataList = entry.getValue();
            final String chromosome = entry.getKey();

            if (geneDataList == null || geneDataList.isEmpty())
                continue;

            LOGGER.debug("calculating fragment size for chromosome({}), fragCount({})", chromosome, mProcessedFragments);

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
                    LOGGER.debug("chromosome({}) processed {} genes, lastGenePos({}) fragCount({}) totalReads({})",
                            chromosome, i, lastGeneEnd, mProcessedFragments, mTotalReadCount);
                }

                lastGeneEnd = geneData.GeneEnd;

                mCurrentReadCount = 0;
                mCurrentGeneData = geneData;
                mRnaBamReader.readBamCounts(GenomeRegions.create(chromosome, geneData.GeneStart, geneData.GeneEnd), this::processBamRead);

                if(mConfig.FragmentLengthsByGene)
                {
                    writeFragmentLengths(geneData);
                    mFragmentLengths.clear();
                }

                if (mProcessedFragments >= mConfig.FragmentLengthMinCount)
                {
                    LOGGER.debug("max fragment length samples reached: {}", mProcessedFragments);
                    break;
                }
            }

            if (mProcessedFragments >= mConfig.FragmentLengthMinCount)
                break;
        }

        if (mConfig.WriteFragmentLengths && !mConfig.FragmentLengthsByGene)
        {
            writeFragmentLengths(null);
        }

        closeBufferedWriter(mWriter);

        if(!mConfig.FragmentLengthsByGene)
        {
            calcSummaryData();
        }
    }

    private void calcSummaryData()
    {
        // work out median or Xth percential etc

        int totalLengths = mFragmentLengths.stream().mapToInt(x -> x[FL_FREQUENCY]).sum();
        int medianCount = totalLengths / 2;
        int percentile75Count = (int) round(totalLengths * 0.75);

        int countsTotal = 0;
        int medianLength = -1;
        int percentile75Length = -1;

        for (final int[] flData : mFragmentLengths)
        {
            if (medianLength == -1 && flData[FL_FREQUENCY] + countsTotal > medianCount)
            {
                medianLength = flData[FL_LENGTH];
            }
            else if (flData[FL_FREQUENCY] + countsTotal > percentile75Count)
            {
                percentile75Length = flData[FL_LENGTH];
                break;
            }

            countsTotal += flData[FL_FREQUENCY];
        }

        LOGGER.debug("fragment length median({}) 75th percentile({})", medianLength, percentile75Length);
    }

    private void processBamRead(@NotNull final SAMRecord record)
    {
        ++mCurrentReadCount;
        ++mTotalReadCount;

        if(mTotalReadCount > 0 && (mTotalReadCount % 10000) == 0)
        {
            LOGGER.trace("currentGene({}:{}) totalReads({})", mCurrentGeneData.GeneId, mCurrentGeneData.GeneName, mTotalReadCount);
        }

        if(mCurrentReadCount >= MAX_GENE_READ_COUNT)
            return;

        if(!isCandidateRecord(record))
            return;

        long posStart = record.getStart();
        long posEnd = record.getEnd();

        for(final TranscriptData transData : mCurrentTransDataList)
        {
            if(transData.exons().stream().anyMatch(x -> !(posStart > x.ExonEnd || posEnd < x.ExonStart)))
                return;
        }

        addFragmentLength(record);
    }

    public void recordFragmentLength(final SAMRecord record, final GeneReadData geneReadData)
    {
        if(!mConfig.WriteFragmentLengths || mConfig.FragmentLengthMinCount > 0)
            return;

        if(!isCandidateRecord(record))
            return;

        long posStart = record.getStart();
        long posEnd = record.getEnd();

        // no part of the read can overlap with any exon
        if(geneReadData.getExonRegions().stream().anyMatch(x -> !(posStart > x.end() || posEnd < x.start())))
            return;

        addFragmentLength(record);
    }

    private boolean isCandidateRecord(final SAMRecord record)
    {
        if(!record.getFirstOfPairFlag() || mRnaBamReader.checkDuplicates(record))
            return false;

        // ignore translocations and inversions
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
            return false;

        // ignore split or unmapped reads
        if(record.getCigar() == null || record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return false;

        return true;
    }

    private void addFragmentLength(final SAMRecord record)
    {
        // long mateStartPos = record.getMateAlignmentStart();

        // long inferredFragLength = mateStartPos > record.getStart() ?
        //        mateStartPos - record.getStart() + record.getReadLength() : record.getEnd() - mateStartPos;

        /* don't discard long reads
        if(inferredFragLength > FRAG_LENGTH_LIMIT)
        {
            // more likely a split exonic read pair
            return;
        }
        */

        int fragmentLength = min(abs(record.getInferredInsertSize()), FRAG_LENGTH_CAP);

        if(fragmentLength == 0)
            return;

        int index = 0;
        boolean exists = false;
        while(index < mFragmentLengths.size())
        {
            final int[] fragLengthCount = mFragmentLengths.get(index);

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
            mFragmentLengths.add(index, newFragLengthCount);
        }

        ++mProcessedFragments;
    }

    public void writeFragmentLengths(@Nullable final EnsemblGeneData geneData)
    {
        if(mConfig.OutputDir.isEmpty() || mFragmentLengths.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("fragment_lengths.csv");
                mWriter = createBufferedWriter(outputFileName, false);

                if(geneData != null)
                {
                    mWriter.write("GeneId,GeneName,Chromosome,");
                }

                mWriter.write("FragmentLength,Count");
                mWriter.newLine();
            }

            for (final int[] fragLengthCount : mFragmentLengths)
            {
                if(geneData != null)
                {
                    mWriter.write(String.format("%s,%s,%s,",
                            geneData.GeneId, geneData.GeneName, geneData.Chromosome));
                }

                mWriter.write(String.format("%d,%d", fragLengthCount[FL_LENGTH], fragLengthCount[FL_FREQUENCY]));
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write fragment length file: {}", e.toString());
        }

    }

}
