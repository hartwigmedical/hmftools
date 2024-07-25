package com.hartwig.hmftools.orange.algo.gripss;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public final class GripssDataLoader
{
    private static final String GRIPSS_FILTERED_SOMATIC_SUFFIX = ".gripss.filtered.somatic.vcf";

    @NotNull
    public static GripssData load(@NotNull String tumorSample, @NotNull String gripssSomaticDir) throws IOException
    {
        String gripsSomaticFilename = vcfPath(tumorSample, gripssSomaticDir);

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());
        VariantContextFilter hotspotFilter = record -> record.hasAttribute("HOTSPOT");
        filter.add(hotspotFilter);

        try
        {
            List<StructuralVariant> sv = StructuralVariantFileLoader.fromFile(gripsSomaticFilename, filter);
            return ImmutableGripssData.builder().allSomaticStructuralVariants(sv).build();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Unable to read structural variant file", e);
        }
    }

    @NotNull
    private static String vcfPath(@NotNull String tumorSample, @NotNull String gripssSomaticDir)
    {
        String filename = gripssSomaticDir + "/" + tumorSample + GRIPSS_FILTERED_SOMATIC_SUFFIX;
        // TODO how does gz extension fallback work in the common framework?
        if(new File(filename).exists())
        {
            return filename;
        }
        return filename + ".gz";
    }
}
