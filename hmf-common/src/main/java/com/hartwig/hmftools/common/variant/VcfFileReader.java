package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class VcfFileReader implements AutoCloseable
{
    private final String mFilename;
    private final boolean mFileValid;
    private final AbstractFeatureReader<VariantContext, LineIterator> mReader;

    private final Map<String, Integer> mGenotypeOrdinals;
    private int mReferenceOrdinal;
    private int mTumorOrdinal;

    private static final Logger LOGGER = LogManager.getLogger(VcfFileReader.class);

    private static final int NO_GENOTYPE_INDEX = -1;

    public VcfFileReader(final String filename)
    {
        this(filename, false);
    }

    public VcfFileReader(final String filename, final boolean requireIndex)
    {
        mFilename = filename;
        mGenotypeOrdinals = Maps.newHashMap();
        mReferenceOrdinal = NO_GENOTYPE_INDEX;
        mTumorOrdinal = NO_GENOTYPE_INDEX;

        if(Files.exists(Paths.get(filename)))
        {
            mReader = getFeatureReader(filename, new VCFCodec(), requireIndex);
            mFileValid = true;

            List<String> vcfSampleNames = ((VCFHeader) mReader.getHeader()).getGenotypeSamples();
            for(int i = 0; i < vcfSampleNames.size(); ++i)
            {
                mGenotypeOrdinals.put(vcfSampleNames.get(i), i);
            }
        }
        else
        {
            mFileValid = false;
            mReader = null;
        }
    }

    public void registerSampleNames(final String referenceId, final String tumorId)
    {
        List<String> vcfSampleNames = ((VCFHeader) mReader.getHeader()).getGenotypeSamples();
        for(int i = 0; i < vcfSampleNames.size(); ++i)
        {
            if(vcfSampleNames.get(i).equals(referenceId))
                mReferenceOrdinal = i;
            else if(vcfSampleNames.get(i).equals(tumorId))
                mTumorOrdinal = i;
        }
    }

    public Map<String, Integer> genotypeOrdinals() { return mGenotypeOrdinals; }
    public int tumorOrdinal() { return mTumorOrdinal; }
    public int referenceOrdinal() { return mReferenceOrdinal; }

    public boolean fileValid() { return mFileValid; }

    public void close()
    {
        try
        {
            if(mReader != null)
                mReader.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to close VCF({}): {}", mFilename, e.toString());
        }
    }

    @Nullable
    public VCFHeader vcfHeader() { return (VCFHeader) mReader.getHeader(); }

    @Nullable
    public AbstractFeatureReader<VariantContext, LineIterator> reader() { return mReader; }

    public CloseableTribbleIterator<VariantContext> iterator()
    {
        try
        {
            return mReader.iterator();
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read variant from file({}): {}", mFilename, e.toString());
            return null;
        }
    }

    @Nullable
    public VariantContext nextVariant() { return iterator().next(); }

    public List<VariantContext> findVariants(final String chromosome, final int lowerBound, final int upperBound)
    {
        CloseableTribbleIterator<VariantContext> iterator = regionIterator(chromosome, lowerBound, upperBound);
        return iterator != null ? iterator.stream().collect(Collectors.toList()) : Collections.emptyList();
    }

    public CloseableTribbleIterator<VariantContext> regionIterator(final ChrBaseRegion region)
    {
        return regionIterator(region.Chromosome, region.start(), region.end());
    }

    public CloseableTribbleIterator<VariantContext> regionIterator(final String chromosome, final int lowerBound, final int upperBound)
    {
        try
        {
            return mReader.query(chromosome, lowerBound, upperBound);
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read variants from file({}): {}", mFilename, e.toString());
            return null;
        }
    }
}
