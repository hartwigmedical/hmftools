package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.TINC_RECOVERED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.TINC_RECOVERED_FLAG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.SageVcfTags.writeTincLevel;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_AVG_READS;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_MAX;
import static com.hartwig.hmftools.sage.tinc.TincAnalyser.filterOutVariant;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.common.variant.pon.PonCache;
import com.hartwig.hmftools.sage.filter.SoftFilter;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class TincVcfWriter
{
    public final TincConfig mConfig;

    public TincVcfWriter(final TincConfig config)
    {
        mConfig = config;
    }

    public static void writeVcf(
            final IndexedFastaSequenceFile refGenome, final String inputVcf, final String outputVcf,
            final List<VariantData> allVariants, final double tincLevel)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(inputVcf);

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        writeTincLevel(vcfHeader, tincLevel);

        if(!vcfHeader.hasFormatLine(TINC_RECOVERED_FLAG))
        {
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(TINC_RECOVERED_FLAG, 0, VCFHeaderLineType.Flag, TINC_RECOVERED_DESC));
        }

        PonCache.addAnnotationHeader(vcfHeader);
        GnomadCache.addAnnotationHeader(vcfHeader);

        VariantVCF vcfFile = new VariantVCF(refGenome, Collections.emptyList(), vcfHeader, outputVcf);

        for(VariantData variant : allVariants)
        {
            if(!variant.isPassing() && !variant.recovered())
            {
                if(filterOutVariant(variant))
                    continue;
            }

            if(variant.gnomadFrequency() != null)
                variant.Context.getCommonInfo().putAttribute(GNOMAD_FREQ, variant.gnomadFrequency());

            if(variant.ponSampleCount() > 0 || variant.ponMaxReadCount() > 0)
            {
                variant.Context.getCommonInfo().putAttribute(PON_COUNT, variant.ponSampleCount());
                variant.Context.getCommonInfo().putAttribute(PON_MAX, variant.ponMaxReadCount());
                variant.Context.getCommonInfo().putAttribute(PON_AVG_READS, variant.ponMeanReadCount());
            }

            if(variant.newFilters() != null && !filtersMatch(variant.newFilters(), variant.Context.getFilters()))
            {
                VariantContext newContext = new VariantContextBuilder(variant.Context)
                        .filters(Collections.emptySet())
                        .make(true);

                if(variant.recovered())
                {
                    newContext.getCommonInfo().putAttribute(TINC_RECOVERED_FLAG, true);
                    newContext.getCommonInfo().addFilter(PASS);
                }
                else
                {
                    variant.newFilters().forEach(x -> newContext.getCommonInfo().addFilter(x.filterName()));
                }

                vcfFile.write(newContext);
            }
            else
            {
                vcfFile.write(variant.Context);
            }
        }

        vcfFile.close();
    }

    public static boolean filtersMatch(final Set<SoftFilter> filters, final Set<String> filterStrings)
    {
        if(filters.size() != filterStrings.size())
            return false;

        return filters.stream().allMatch(x -> filterStrings.contains(x.filterName()));
    }

}
