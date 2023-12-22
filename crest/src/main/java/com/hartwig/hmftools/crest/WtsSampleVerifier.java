package com.hartwig.hmftools.crest;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

// Proof of concept - check if wts sample is from same patient
// as wgs sample.
public class WtsSampleVerifier
{

    public static void main(String[] args) throws IOException
    {
        String rnaAnnotatedGermlineVcf = args[0];
        String rnaSample = args[1];
        check(rnaAnnotatedGermlineVcf, rnaSample);
    }

    public static double check(String germlineVcf, String rnaSample) throws IOException
    {
        int supported = 0;
        var total = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(germlineVcf, new VCFCodec(), false))
        {
            for(VariantContext context : reader.iterator())
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(decorator.filter().equals("PASS") && decorator.type() == VariantType.SNP && !decorator.gene().isEmpty())
                {
                    AllelicDepth rnaDepth = decorator.allelicDepth(rnaSample);
                    if(rnaDepth.totalReadCount() >= 10)
                    {
                        total += 1;
                        if(rnaDepth.alleleReadCount() >= 1)
                        {
                            supported += 1;
                        }
                    }
                }
            }
        }
        double fraction = total > 0 ? supported * 1D / total : 0D;
        System.out.println("Supported: " + supported + " Total: " + total + " Fraction: " + fraction);
        return fraction;
    }
}