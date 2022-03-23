package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.purple.somatic.VariantHotspotEnrichment;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.purple.StructuralVariantCache;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFFileReader;

public class SampleData
{
    public final String ReferenceId;
    public final String SampleId;

    public final AmberData Amber;
    public final CobaltData Cobalt;
    public final StructuralVariantCache SvCache;
    public final List<SomaticVariant> FittingSomaticVariants;

    public SampleData(final String referenceId, final String sampleId,
            final AmberData amber, final CobaltData cobalt, final StructuralVariantCache svCache)
    {
        ReferenceId = referenceId;
        SampleId = sampleId;

        Amber = amber;
        Cobalt = cobalt;
        SvCache = svCache;
        FittingSomaticVariants = Lists.newArrayList();
    }

    public void loadSomatics(final String somaticVcf, final ReferenceData referenceData)
    {
        if(somaticVcf.isEmpty())
            return;

        final SomaticVariantFactory factory = new SomaticVariantFactory(new PassingVariantFilter(), new SGTFilter());

        final VariantHotspotEnrichment hotspotEnrichment = new VariantHotspotEnrichment(referenceData.SomaticHotspots, true);

        PPL_LOGGER.info("loading somatic variants from {}", somaticVcf);

        try (VCFFileReader vcfReader = new VCFFileReader(new File(somaticVcf), false))
        {
            for(VariantContext variantContext : vcfReader)
            {
                if(factory.test(variantContext))
                {
                    Optional<SomaticVariant> variant = factory.createVariant(SampleId, hotspotEnrichment.processVariant(variantContext));

                    if(variant.isPresent() && variant.get().isSnp())
                        FittingSomaticVariants.add(variant.get());
                }
            }
        }
    }
}
