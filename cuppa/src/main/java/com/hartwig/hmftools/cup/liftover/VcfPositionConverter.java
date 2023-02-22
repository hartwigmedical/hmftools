package com.hartwig.hmftools.cup.liftover;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.isKnownIndel;
import static com.hartwig.hmftools.cup.liftover.SnvLiftover.LIFTOVER_FILE;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.Callable;

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
    private final CoordMappingCache mMappingCache;

    private BufferedWriter mWriter;
    private int mMappingIndex;
    private final boolean mMappingEnabled;

    private static final int UNMAPPED_POSITION = -1;

    public VcfPositionConverter(
            final String sampleId, final String vcfFile, final CoordMappingCache mappingCache, final LiftoverConfig config)
    {
        mSampleId = sampleId;
        mConfig = config;

        mVcfFile = vcfFile;
        mMappingCache = mappingCache;
        mMappingEnabled = mMappingCache.hasMappings();
        mMappingIndex = 0;

        mOutputFile = mConfig.OutputDir + mSampleId + LIFTOVER_FILE;
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

                if(mConfig.ApplyFilters)
                {
                    if(variant.Type == VariantType.MNP)
                        continue;

                    if(variant.Type == VariantType.INDEL && !isKnownIndel(variant.Gene, variant.RepeatCount, variant.Type))
                        continue;
                }

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

    private int convertPosition(final String chromosome, final int position)
    {
        for(; mMappingIndex < mMappingCache.getMappings().size(); ++mMappingIndex)
        {
            CoordMapping mapping = mMappingCache.getMappings().get(mMappingIndex);

            if(!mapping.Chromosome.equals(chromosome))
                continue;

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
            writer.write("Chromosome,Position,Ref,Alt,Type,RepeatCount,Gene,TriNucContext,PrevPosition");
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

            mWriter.write(format("%s,%d,%s,%s,%s,%d,%s,%s,%d",
                    outputChr, newPosition, variant.Ref, variant.Alt, variant.Type,
                    variant.RepeatCount, variant.Gene, variant.TrinucleotideContext, variant.Position));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to write variant: {}", e.toString());
        }
    }
}
