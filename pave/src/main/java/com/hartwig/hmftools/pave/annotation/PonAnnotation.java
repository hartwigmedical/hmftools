package com.hartwig.hmftools.pave.annotation;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.pon.PonCache;
import com.hartwig.hmftools.common.variant.pon.PonChrCache;
import com.hartwig.hmftools.common.variant.pon.PonVariantData;
import com.hartwig.hmftools.pave.VariantData;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PonAnnotation extends AnnotationData implements Callable
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
    public Long call()
    {
        for(String chromosome : mInitialChromosomes)
        {
            mPonCache.loadPonEntries(chromosome);
        }

        return (long)0;
    }

    public void annotateVariant(final VariantData variant, final PonChrCache chrCache)
    {
        PonVariantData ponData = chrCache.getPonData(variant.Position, variant.Ref, variant.Alt);
        if(ponData == null)
            return;

        variant.setPonFrequency(ponData.Samples, ponData.MaxSampleReads, ponData.meanReadCount());
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
