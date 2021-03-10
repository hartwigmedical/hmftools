package com.hartwig.hmftools.common.drivercatalog.panel;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

final class GermlineResources {

    static List<VariantContext> grch37Whitelist() throws IOException {
        return getResource(getResourceURL("/drivercatalog/GermlineHotspots.whitelist.hg19.vcf"));
    }

    static List<VariantContext> grch38Whitelist() throws IOException {
        return getResource(getResourceURL("/drivercatalog/GermlineHotspots.whitelist.hg38.vcf"));
    }

    static List<VariantContext> grch37Blacklist() throws IOException {
        return getResource(getResourceURL("/drivercatalog/GermlineHotspots.blacklist.hg19.vcf"));
    }

    static List<VariantContext> grch38Blacklist() throws IOException {
        return getResource(getResourceURL("/drivercatalog/GermlineHotspots.blacklist.hg38.vcf"));
    }

    private static List<VariantContext> getResource(@NotNull final String file) throws IOException {
        return getFeatureReader(file, new VCFCodec(), false).iterator().toList();
    }

    private static String getResourceURL(@NotNull String location) {
        return DriverGenePanelConversion.class.getResource(location).toString();
    }
}
