package com.hartwig.hmftools.geneutils.drivers;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class GermlineResources
{
    static List<VariantContext> whitelist37() throws IOException
    {
        return resource(resourceURL("/drivers/GermlineHotspots.whitelist.37.vcf"));
    }

    static List<VariantContext> whitelist38() throws IOException
    {
        return resource(resourceURL("/drivers/GermlineHotspots.whitelist.38.vcf"));
    }

    static List<VariantContext> blacklist37() throws IOException
    {
        return resource(resourceURL("/drivers/GermlineHotspots.blacklist.37.vcf"));
    }

    static List<VariantContext> blacklist38() throws IOException
    {
        return resource(resourceURL("/drivers/GermlineHotspots.blacklist.38.vcf"));
    }

    private static List<VariantContext> resource(final String file) throws IOException
    {
        return getFeatureReader(file, new VCFCodec(), false).iterator().toList();
    }

    private static String resourceURL(final String location)
    {
        return GenerateDriverGeneFiles.class.getResource(location).toString();
    }
}
