package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_AVG_READS;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_MAX;

import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.pon.PonCache;
import com.hartwig.hmftools.common.variant.pon.PonChrCache;
import com.hartwig.hmftools.common.variant.pon.PonVariantData;
import com.hartwig.hmftools.pave.VariantData;

public class PonAnnotation extends AnnotationData implements Callable<Void>
{
    private final String mPonFilename;

    private final PonCache mPonCache;

    public static final String PON_ARTEFACT_FILTER = "PONArtefact";

    public PonAnnotation(final String filename, boolean loadOnDemand)
    {
        mPonFilename = filename;
        mPonCache = new PonCache(filename, loadOnDemand);
    }

    @Override
    public String type() { return "PON"; }

    @Override
    public boolean enabled() { return mPonFilename != null; }

    @Override
    public boolean hasValidData() { return mPonCache.hasValidData(); }

    public boolean loadFilters(final String filtersConfig)
    {
        if(filtersConfig == null)
            return true;

        return mPonCache.loadFilters(filtersConfig);
    }

    public synchronized PonChrCache getChromosomeCache(final String chromosome)
    {
        return mPonCache.getChromosomeCache(chromosome);
    }

    @Override
    public synchronized void onChromosomeComplete(final String chromosome)
    {
        mPonCache.removeCompleteChromosome(chromosome);
    }

    @Override
    public Void call()
    {
        for(String chromosome : mInitialChromosomes)
        {
            mPonCache.loadPonEntries(chromosome);
        }

        return null;
    }

    public void annotateVariant(final VariantData variant, final PonChrCache chrCache)
    {
        PonVariantData ponData = chrCache.getPonData(variant.Position, variant.Ref, variant.Alt);
        if(ponData == null)
            return;

        variant.setPonFrequency(ponData.Samples, ponData.MaxSampleReads, ponData.meanReadCount());
    }

    public static void annotateFromContext(final VariantData variant)
    {
        if(variant.context().hasAttribute(PON_COUNT))
        {
            variant.setPonFrequency(
                    variant.context().getAttributeAsInt(PON_COUNT, 0),
                    variant.context().getAttributeAsInt(PON_MAX, 0),
                    variant.context().getAttributeAsInt(PON_AVG_READS, 0));
        }
    }

    public boolean filterOnTierCriteria(final VariantTier tier, final int ponSampleCount, final int ponMaxSampleReads)
    {
        return mPonCache.filterOnTierCriteria(tier, ponSampleCount, ponMaxSampleReads);
    }

    public boolean hasEntry(final String chromosome, final int position, final String ref, final String alt)
    {
        return mPonCache.hasEntry(chromosome, position, ref, alt);
    }
}
