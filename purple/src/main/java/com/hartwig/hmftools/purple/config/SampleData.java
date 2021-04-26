package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.File;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.enrich.VariantHotspotEnrichment;
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
    public final List<SomaticVariant> SomaticVariants;
    public final List<SomaticVariant> FittingSomaticVariants;

    public SampleData(final String referenceId, final String sampleId,
            final AmberData amber, final CobaltData cobalt, final StructuralVariantCache svCache)
    {
        ReferenceId = referenceId;
        SampleId = sampleId;

        Amber = amber;
        Cobalt = cobalt;
        SvCache = svCache;
        SomaticVariants = Lists.newArrayList();
        FittingSomaticVariants = Lists.newArrayList();
    }

    public void loadSomatics(final String somaticVcf, final ReferenceData referenceData, boolean tumorOnlyMode)
    {
        if(somaticVcf.isEmpty())
        {
            PPL_LOGGER.info("Somatic variants support disabled.");
            return;
        }

        final SomaticVariantFactory factory = new SomaticVariantFactory(new PassingVariantFilter(), new SGTFilter());
        final Consumer<VariantContext> convertThenSave = context -> factory.createVariant(SampleId, context).ifPresent(SomaticVariants::add);

        final VariantHotspotEnrichment enrich = new VariantHotspotEnrichment(referenceData.SomaticHotspots, convertThenSave);

        PPL_LOGGER.info("Loading somatic variants from {}", somaticVcf);

        try (VCFFileReader vcfReader = new VCFFileReader(new File(somaticVcf), false))
        {
            for(VariantContext variantContext : vcfReader)
            {
                if(factory.test(variantContext))
                {
                    enrich.accept(variantContext);
                }
            }
        }

        enrich.flush();

        if(!tumorOnlyMode)
        {
            SomaticVariants.stream().filter(SomaticVariant::isSnp).forEach(x -> FittingSomaticVariants.add(x));
        }
    }


}
