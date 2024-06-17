package com.hartwig.hmftools.cup.liftover;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.CoordMapping;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.cup.somatics.SomaticVariant;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class VcfPositionConverter implements Callable
{
    private final LiftoverConfig mConfig;
    private final String mSampleId;
    private final String mVcfFile;
    private final String mOutputFile;
    private final GenomeLiftoverCache mMappingCache;
    private final boolean mMappingEnabled;

    private BufferedWriter mWriter;

    // look-up state
    private int mCurentMappingIndex;
    private String mCurrentMappingChromosome;
    private List<CoordMapping> mChromosomeMappings;

    public VcfPositionConverter(
            final String sampleId, final String vcfFile, final LiftoverConfig config)
    {
        mSampleId = sampleId;
        mConfig = config;

        mVcfFile = vcfFile;
        mMappingCache = new GenomeLiftoverCache(true);
        mMappingEnabled = mMappingCache.hasMappings();
        mCurentMappingIndex = 0;
        mCurrentMappingChromosome = "";
        mChromosomeMappings = null;

        mOutputFile = SomaticVariant.generateFilename(mConfig.OutputDir, mSampleId);
        mWriter = null;
    }

    @Override
    public Long call()
    {
        if(mConfig.KeepExisting && Files.exists(Paths.get(mOutputFile)))
        {
            CUP_LOGGER.info("sample({}) output exists, skipping", mSampleId);
            return (long)0;
        }

        mWriter = initialiseWriter();

        try
        {
            int variantCount = 0;

            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(mVcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                if(mConfig.ApplyFilters && variantContext.isFiltered())
                    continue;

                SomaticVariant variant = SomaticVariant.fromContext(variantContext);

                if(mConfig.ApplyFilters && variant.Type == VariantType.MNP)
                    continue;

                int convertedPosition = mMappingEnabled ? convertPosition(variant.Chromosome, variant.Position) : variant.Position;
                writeVariant(variant, convertedPosition);
                ++variantCount;
            }

            CUP_LOGGER.debug("sample({}) converted {} variants", mSampleId, variantCount);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to read somatic VCF file({}): {}", mSampleId, mVcfFile, e.toString());
        }

        closeBufferedWriter(mWriter);

        return (long)0;
    }

    @VisibleForTesting
    public int convertPosition(final String chromosome, final int position)
    {
        if(!mCurrentMappingChromosome.equals(chromosome))
        {
            mChromosomeMappings = mMappingCache.getChromosomeMappings(chromosome, true);
            mCurrentMappingChromosome = chromosome;
            mCurentMappingIndex = 0;
        }

        for(; mCurentMappingIndex < mChromosomeMappings.size(); ++mCurentMappingIndex)
        {
            CoordMapping mapping = mChromosomeMappings.get(mCurentMappingIndex);

            if(mapping.SourceEnd < position)
                continue;

            if(mapping.SourceStart > position)
                return UNMAPPED_POSITION;

            return mapping.convertPosition(position);
        }

        return UNMAPPED_POSITION;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);
            writer.write(format("%s,PrevPosition", SomaticVariant.csvHeader()));
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to create output file: {}", e.toString());
            return null;
        }
    }

    private void writeVariant(final SomaticVariant variant, int newPosition)
    {
        if(mWriter == null)
            return;

        try
        {
            String outputChr = mMappingEnabled ? RefGenomeVersion.V38.versionedChromosome(variant.Chromosome) : variant.Chromosome;

            mWriter.write(format("%s,%d,%s,%s,%s,%s,%s,%d,%d",
                    outputChr, newPosition, variant.Ref, variant.Alt, variant.Type,
                    variant.Gene, variant.TrinucleotideContext, variant.RepeatCount, variant.Position));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to write variant: {}", e.toString());
        }
    }
}
