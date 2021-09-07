package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svtools.pon.PonBuilder.PON_LOGGER;

import java.io.IOException;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class PonSampleTask implements Callable
{
    private final int mTaskId;
    private final Map<String,String> mSampleVcfFiles;
    private final PonStore mPonStore;

    private StructuralVariantFactory mSvFactory;
    private int mProcessedVariants;

    public PonSampleTask(final int taskId, final PonStore ponStore)
    {
        mTaskId = taskId;
        mPonStore = ponStore;

        mSampleVcfFiles = Maps.newHashMap();

        mSvFactory = null;
        mProcessedVariants = 0;
    }

    @Override
    public Long call()
    {
        processSamples();
        return (long)0;
    }

    public Map<String,String> getSampleVcfFiles() { return mSampleVcfFiles; }

    public void processSamples()
    {
        for(Map.Entry<String,String> entry : mSampleVcfFiles.entrySet())
        {
            processVcf(entry.getKey(), entry.getValue());
        }
    }

    private void processVcf(final String sampleId, final String vcfFile)
    {
        clearSampleData();

        try
        {
            PON_LOGGER.info("{}: processing sampleId({}) vcfFile({})", mTaskId, sampleId, vcfFile);

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));

        }
        catch(IOException e)
        {
            PON_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }

        PON_LOGGER.info("{}: sample({}) read {} variants", mTaskId, sampleId, mProcessedVariants);

        // log current PON stats
        PON_LOGGER.info("{}: {}}", mTaskId, mPonStore.statsString());
    }

    private void clearSampleData()
    {
        mProcessedVariants = 0;
        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
    }

    private void processVariant(final VariantContext variant)
    {
        PON_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            // PON_LOGGER.debug("sample({}) processed {} variants, VCF-unmatched({})",
            //        mConfig.SampleId, mProcessedVariants, mSvFactory.unmatched().size());
        }

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // wait for both breakends to be added
        if(currentSvCount == mSvFactory.results().size())
            return;

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        if(sv.type() == SGL)
        {
            mPonStore.addLocation(sv.chromosome(true), sv.orientation(true), sv.position(true).intValue());
        }
        else
        {
            mPonStore.addLocation(
                    sv.chromosome(true), sv.chromosome(false),
                    sv.orientation(true), sv.orientation(false),
                    sv.position(true).intValue(), sv.position(false).intValue());
        }
    }

    private final StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

}
