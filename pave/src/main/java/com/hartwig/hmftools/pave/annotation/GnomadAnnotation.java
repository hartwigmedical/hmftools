package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.variant.pon.GnomadCache.GNOMAD_FREQUENCY_DIR;
import static com.hartwig.hmftools.common.variant.pon.GnomadCache.GNOMAD_FREQUENCY_FILE;
import static com.hartwig.hmftools.common.variant.pon.GnomadCache.GNOMAD_NO_FILTER;

import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.pon.GnomadChrCache;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.pave.VariantData;

public class GnomadAnnotation extends AnnotationData implements Callable
{
    private final GnomadCache mGnomadCache;
    private final boolean mNoFilter;

    public GnomadAnnotation(final ConfigBuilder configBuilder)
    {
        mGnomadCache = new GnomadCache(
                RefGenomeVersion.from(configBuilder),
                configBuilder.getValue(GNOMAD_FREQUENCY_FILE),
                configBuilder.getValue(GNOMAD_FREQUENCY_DIR));

        mNoFilter = configBuilder.hasFlag(GNOMAD_NO_FILTER);
    }

    @Override
    public String type() { return "Gnomad frequency"; }

    @Override
    public boolean enabled() { return mGnomadCache.enabled(); }

    public boolean applyFilter() { return !mNoFilter; }

    @Override
    public boolean hasValidData() { return mGnomadCache.hasValidData(); }

    public void annotateVariant(final VariantData variant, final GnomadChrCache chrCache)
    {
        Double gnomadFreq = chrCache.getFrequency(variant.isMnv(), variant.Ref, variant.Alt, variant.Position);

        if(gnomadFreq != null)
        {
            variant.setGnomadFrequency(gnomadFreq);
        }
    }

    public GnomadChrCache getChromosomeCache(final String chromosome)
    {
        return mGnomadCache.getChromosomeCache(chromosome);
    }

    @Override
    public void onChromosomeComplete(final String chromosome)
    {
        mGnomadCache.removeCompleteChromosome(chromosome);
    }

    @Override
    public Long call()
    {
        mGnomadCache.initialise(mInitialChromosomes);

        return (long)0;
    }
}
