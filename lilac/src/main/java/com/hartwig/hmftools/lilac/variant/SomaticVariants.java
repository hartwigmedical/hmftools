package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.List;
import java.util.Set;

public class SomaticVariants
{
    private final Set<CodingEffect> UNKNOWN_CODING_EFFECT;
    private final LilacConfig mConfig;

    public SomaticVariants(final LilacConfig config)
    {
        mConfig = config;
        UNKNOWN_CODING_EFFECT = Sets.newHashSet(CodingEffect.NONE, CodingEffect.UNDEFINED);
    }

    public final List<VariantContextDecorator> readSomaticVariants()
    {
        if(mConfig.TumorBam.isEmpty() || mConfig.SomaticVcf.isEmpty())
            return Lists.newArrayList();

        final List<VariantContextDecorator> results = Lists.newArrayList();

        String chromosome = HLA_TRANSCRIPTS.get(0).chromosome();
        int minPosition = HLA_TRANSCRIPTS.stream().mapToInt(x -> (int)x.start()).min().orElse(0);
        int maxPosition = HLA_TRANSCRIPTS.stream().mapToInt(x -> (int)x.start()).max().orElse(0);

        LL_LOGGER.info("Reading somatic vcf: ", mConfig.SomaticVcf);

        VCFFileReader fileReader = new VCFFileReader(new File(mConfig.SomaticVcf), false);

        final CloseableIterator<VariantContext> variantIter = fileReader.isQueryable() ?
                fileReader.query(chromosome, minPosition, maxPosition) : fileReader.iterator();

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

        LL_LOGGER.info("    found ${result.size} HLA variants", 0);
        return results;

        /*
        this.logger.info("Reading somatic vcf: " + this.mConfig.getSomaticVcf());
        List result = new ArrayList();
        VCFFileReader fileReader = new VCFFileReader(new File(this.mConfig.getSomaticVcf()), false);
        CloseableIterator closeableIterator = iterator =
                fileReader.isQueryable() ? fileReader.query(contig, (int) minPosition, (int) maxPosition) : fileReader.iterator();
        Intrinsics.checkExpressionValueIsNotNull((Object) closeableIterator, (String) "iterator");
        Iterator iterator2 = it = (Iterator) closeableIterator;
        while(iterator2.hasNext())
        {
            VariantContext variantContext = (VariantContext) iterator2.next();
            VariantContextDecorator enriched = new VariantContextDecorator(variantContext);
            if(!LilacApplication.Companion.getHLA_GENES().contains(enriched.gene())
                    || this.UNKNOWN_CODING_EFFECT.contains(enriched.canonicalCodingEffect()) || !enriched.isPass())
            {
                continue;
            }
            result.add(enriched);
        }
        fileReader.close();
        this.logger.info("    found " + result.size() + " HLA variants");
        return result;

         */
    }

}
