package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Germline SNV or INDEL.
public class GermlineMutation extends Variant
{
    private final GermlineVariant mVariant;

    private GermlineMutation(final GermlineVariant variant)
    {
        mVariant = variant;
    }

    private static final Logger LOGGER = LogManager.getLogger(GermlineMutation.class);

    // TODO: only select if driver = true
    @Override
    public boolean isDriver()
    {
        return mVariant.reported();
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        return buildMutationProbe(mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(), PROBE_LENGTH, refGenome);
    }

    @Override
    public List<ProximateLocations.Location> checkedLocations()
    {
        return List.of(new ProximateLocations.Location(mVariant.chromosome(), mVariant.position()));
    }

    @Override
    public String toString()
    {
        return format("%s %s:%s %s>%s",
                mVariant.type(), mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
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
