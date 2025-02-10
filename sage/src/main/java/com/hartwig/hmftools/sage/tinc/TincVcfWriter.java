package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TINC_LEVEL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TINC_RECOVERED_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TINC_RECOVERED_FLAG;
import static com.hartwig.hmftools.sage.tinc.TincAnalyser.filterOutVariant;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
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

        if(tincLevel > 0)
        {
            vcfHeader.addMetaDataLine(new VCFHeaderLine(TINC_LEVEL, String.valueOf(tincLevel)));
        }

        if(!vcfHeader.hasFormatLine(TINC_RECOVERED_FLAG))
        {
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(TINC_RECOVERED_FLAG, 0, VCFHeaderLineType.Flag, TINC_RECOVERED_DESC));
        }

        VariantVCF vcfFile = new VariantVCF(refGenome, Collections.emptyList(), vcfHeader, outputVcf);

        for(VariantData variant : allVariants)
        {
            if(!variant.isPassing() && !variant.recovered())
            {
                if(filterOutVariant(variant))
                    continue;
            }

            if(variant.newFilters() != null && variant.newFilters().size() != variant.Context.getFilters().size())
            {
                VariantContext newContext = new VariantContextBuilder(variant.Context).make(true);

                if(variant.recovered())
                {
                    newContext.getCommonInfo().putAttribute(TINC_RECOVERED_FLAG, true);
                    newContext.getCommonInfo().addFilter(PASS);
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
}
