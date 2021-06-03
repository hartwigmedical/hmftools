package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class SomaticVariantFinder
{
    private final Set<CodingEffect> UNKNOWN_CODING_EFFECT;
    private final LilacConfig mConfig;
    private final Map<String,TranscriptData> mHlaTranscriptData;

    public SomaticVariantFinder(final LilacConfig config, final Map<String, TranscriptData> transcriptData)
    {
        mConfig = config;
        mHlaTranscriptData = transcriptData;
        UNKNOWN_CODING_EFFECT = Sets.newHashSet(CodingEffect.NONE, CodingEffect.UNDEFINED);
    }

    public List<VariantContextDecorator> readSomaticVariants()
    {
        if(mConfig.TumorBam.isEmpty() || mConfig.SomaticVariantsFile.isEmpty())
            return Lists.newArrayList();

        final List<VariantContextDecorator> results = Lists.newArrayList();

        int minPosition = mHlaTranscriptData.values().stream().mapToInt(x -> x.TransStart).min().orElse(0);
        int maxPosition = mHlaTranscriptData.values().stream().mapToInt(x -> x.TransEnd).max().orElse(0);

        LL_LOGGER.info("reading somatic vcf: ", mConfig.SomaticVariantsFile);

        VCFFileReader fileReader = new VCFFileReader(new File(mConfig.SomaticVariantsFile), false);

        final CloseableIterator<VariantContext> variantIter = fileReader.isQueryable() ?
                fileReader.query(HLA_CHR, minPosition, maxPosition) : fileReader.iterator();

        while(variantIter.hasNext())
        {
            VariantContext variant = variantIter.next();
            VariantContextDecorator enriched = new VariantContextDecorator(variant);

            if (HLA_GENES.contains(enriched.gene()) && enriched.isPass()
            && !UNKNOWN_CODING_EFFECT.contains(enriched.canonicalCodingEffect()))
            {
                results.add(enriched);
            }
        }

        fileReader.close();

        LL_LOGGER.info("  found {} HLA variants", results.size());
        return results;
    }

}
