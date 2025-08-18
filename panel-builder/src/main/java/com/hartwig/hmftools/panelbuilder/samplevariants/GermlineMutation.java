package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.panelbuilder.ProbeTarget;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Germline SNV or INDEL.
public class GermlineMutation implements Variant
{
    private final GermlineVariant mVariant;

    private GermlineMutation(final GermlineVariant variant)
    {
        mVariant = variant;
    }

    private static final Logger LOGGER = LogManager.getLogger(GermlineMutation.class);

    // TODO: use in extra info
    private String gene()
    {
        return mVariant.gene();
    }

    @Override
    public boolean isDriver()
    {
        return mVariant.reported();
    }

    @Override
    public ProbeTarget generateProbeTarget()
    {
        return buildMutationProbe(mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(), PROBE_LENGTH);
    }

    @Override
    public TargetMetadata.Type targetType()
    {
        if(isDriver())
        {
            return TargetMetadata.Type.SAMPLE_GERMLINE_SNV_INDEL_DRIVER;
        }
        else
        {
            // Shouldn't happen because nondrivers are filtered out.
            throw new IllegalStateException("Unhandled germline mutation type");
        }
    }

    @Override
    public String toString()
    {
        return format("%s:%s %s>%s %s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(), mVariant.type());
    }

    public static List<GermlineMutation> load(final String sampleId, final String purpleDir)
    {
        String vcfFile = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);

        List<GermlineVariant> germlineVariants;
        try
        {
            germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);
        }
        catch(Exception e)
        {
            throw new RuntimeException(format("Failed to load germline variants from file: %s", vcfFile));
        }

        List<GermlineMutation> variants = germlineVariants.stream().map(GermlineMutation::new).toList();

        LOGGER.info("Loaded {} germline mutations", sampleId);
        variants.forEach(variant -> LOGGER.trace("GermlineMutation: {}", variant));

        return variants;
    }
}
