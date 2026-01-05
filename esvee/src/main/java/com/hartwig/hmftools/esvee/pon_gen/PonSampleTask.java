package com.hartwig.hmftools.esvee.pon_gen;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.VariantContext;

public class PonSampleTask implements Callable<Void>
{
    private final List<String> mSampleVcfFiles;
    private final PonConfig mConfig;
    private final PonStore mPonStore;

    private StructuralVariantFactory mSvFactory;
    private int mProcessedVariants;

    private final Set<String> mUniqueSVs;
    private final Set<String> mUniqueSGLs;

    public PonSampleTask(final PonConfig config, final PonStore ponStore)
    {
        mConfig = config;
        mPonStore = ponStore;

        mSampleVcfFiles = Lists.newArrayList();

        mSvFactory = null;
        mProcessedVariants = 0;

        mUniqueSVs = Sets.newHashSet();
        mUniqueSGLs = Sets.newHashSet();
    }

    public void addSampleVcf(final String sampleVcf) { mSampleVcfFiles.add(sampleVcf); }

    @Override
    public Void call()
    {
        for(String sampleVcf : mSampleVcfFiles)
        {
            processVcf(sampleVcf);
        }

        return null;
    }

    private void processVcf(final String sampleVcf)
    {
        clearSampleData();

        VcfFileReader vcfReader = new VcfFileReader(sampleVcf);

        for(VariantContext variantContext : vcfReader.iterator())
        {
            if(mConfig.FilterOnPass && !variantContext.getFilters().isEmpty() && !variantContext.getFilters().contains(PASS_FILTER))
                continue;

            processVariant(variantContext);
        }

        SV_LOGGER.debug("sampleVcf({}) read {} variants", sampleVcf, mProcessedVariants);
    }

    private void clearSampleData()
    {
        mProcessedVariants = 0;
        mSvFactory = StructuralVariantFactory.build(new AlwaysPassFilter());
        mUniqueSVs.clear();
        mUniqueSGLs.clear();
    }

    private static final int VARIANT_LOG_COUNT = 10_000;

    private void processVariant(final VariantContext variant)
    {
        SV_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % VARIANT_LOG_COUNT) == 0)
        {
            // SV_LOGGER.debug("sample({}) processed {} variants, VCF-unmatched({})",
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

        // ignore deduplicated SVs
        if(sv.filter().contains(FilterType.DUPLICATE.vcfTag()))
            return;

        if(!checkAddNew(sv))
            return;

        if(sv.type() == SGL)
        {
            mPonStore.addSgl(sv.chromosome(true), sv.orientation(true), sv.position(true).intValue());
        }
        else
        {
            mPonStore.addSv(
                    sv.chromosome(true), sv.chromosome(false),
                    sv.orientation(true), sv.orientation(false),
                    sv.position(true).intValue(), sv.position(false).intValue());
        }
    }

    private boolean checkAddNew(final StructuralVariant sv)
    {
        if(mConfig.SpecificChrRegions.hasFilters())
        {
            if(mConfig.SpecificChrRegions.excludePosition(sv.chromosome(true), sv.position(true)))
                return false;

            if(sv.type() != SGL && mConfig.SpecificChrRegions.excludePosition(sv.chromosome(false), sv.position(false)))
                return false;
        }

        String id = format("%s_%d_%d", sv.chromosome(true), sv.position(true), sv.orientation(true));

        if(sv.type() == SGL)
        {
            if(mUniqueSGLs.contains(id))
                return false;

            mUniqueSGLs.add(id);
            return true;
        }

        id = format("%s_%s_%d_%d", id, sv.chromosome(false), sv.position(false), sv.orientation(false));

        if(mUniqueSVs.contains(id))
            return false;

        mUniqueSVs.add(id);
        return true;
    }

    private StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }
}
