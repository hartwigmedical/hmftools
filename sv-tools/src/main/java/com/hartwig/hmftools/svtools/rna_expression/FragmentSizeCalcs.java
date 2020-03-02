package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.from;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// calculate fragment length distribution for a sample
// this can be done either independently from fragment counting or lengths can be registering during that process

public class FragmentSizeCalcs
{
    private final RnaExpConfig mConfig;
    private final SvGeneTranscriptCollection mGeneTransCache;

    private final List<int[]> mFragmentLengths;
    private BufferedWriter mWriter;

    private EnsemblGeneData mCurrentGeneData;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentReadCount;
    private int mTotalReadCount;
    private int mProcessedFragments;

    private static final Logger LOGGER = LogManager.getLogger(FragmentSizeCalcs.class);

    public FragmentSizeCalcs(final RnaExpConfig config, final SvGeneTranscriptCollection geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mCurrentGeneData = null;
        mCurrentTransDataList = null;
        mCurrentReadCount = 0;
        mTotalReadCount = 0;
        mProcessedFragments = 0;
        mWriter = writer;

        mFragmentLengths = Lists.newArrayList();
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 20;
    private static final int MAX_TRAN_EXONS = 20;
    private static final int MAX_GENE_READ_COUNT = 1000; // to avoid impact of highly enriched genes

    private static int FRAG_LENGTH_CAP = 3000; // to prevent map blowing out in size

    public void calcSampleFragmentSize()
    {
        // measure fragment lengths for non-split reads in purely intronic regions
        if(mGeneTransCache == null)
        {
            LOGGER.error("fragment length calculator missing gene cache");
            return;
        }

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        BamSlicer bamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        // walk through each chromosome, ignoring any gene which overlaps the previous gene
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
                List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, geneData.GeneStart, geneData.GeneEnd));

                bamSlicer.slice(samReader, regions, this::processBamRead);

                if(mConfig.FragmentLengthsByGene)
                {
                    writeFragmentLengths(mWriter, mFragmentLengths, geneData);
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
            writeFragmentLengths(mWriter, mFragmentLengths, null);
        }
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

    public void addFragmentLength(final SAMRecord record)
    {
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

    public void setConfigLengthDistribution()
    {
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
        }
    }

    public static BufferedWriter createFragmentLengthWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("fragment_lengths.csv");

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
            LOGGER.error("failed to write fragment length file: {}", e.toString());
            return null;
        }
    }

    public void writeFragmentLengths(final EnsemblGeneData geneData)
    {
        writeFragmentLengths(mWriter,mFragmentLengths, mConfig.FragmentLengthsByGene ? geneData : null);
    }

    public synchronized static void writeFragmentLengths(
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
            LOGGER.error("failed to write fragment length file: {}", e.toString());
        }

    }

}
