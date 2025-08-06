package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class GermlineMutation extends Variant
{
    private final GermlineVariant mVariant;

    public GermlineMutation(final GermlineVariant variant)
    {
        mVariant = variant;
    }

    private static final Logger LOGGER = LogManager.getLogger(GermlineMutation.class);

    @Override
    public CategoryType categoryType()
    {
        return CategoryType.GERMLINE_MUTATION;
    }

    @Override
    public String description()
    {
        return format("%s:%s %s>%s %s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(), mVariant.type());
    }

    @Override
    public String gene()
    {
        return mVariant.gene();
    }

    @Override
    public double copyNumber()
    {
        return mVariant.adjustedCopyNumber();
    }

    @Override
    public double vaf()
    {
        return mVariant.adjustedVAF();
    }

    @Override
    public int tumorFragments()
    {
        return mVariant.allelicDepth().AlleleReadCount;
    }

    @Override
    public boolean reported()
    {
        return true;
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        return buildMutationProbe(mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(), PROBE_LENGTH, refGenome);
    }

    @Override
    public boolean checkFilters()
    {
        return false;
    }

    @Override
    public boolean checkAndRegisterLocation(ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariant.chromosome(), mVariant.position()))
        {
            return false;
        }

        registeredLocations.addRegisteredLocation(mVariant.chromosome(), mVariant.position());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<GermlineMutation> load(final String sampleId, final String purpleDir)
    {
        ArrayList<GermlineMutation> variants = new ArrayList<>();

        // TODO: unused - bug?
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);

        try
        {
            List<GermlineVariant> germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);

            germlineVariants.stream()
                    .filter(GermlineVariant::reported)
                    .forEach(v -> variants.add(new GermlineMutation(v)));
        }
        catch(Exception e)
        {
            throw new RuntimeException(format("Failed to load germline variants from file: %s", vcfFile));
        }

        LOGGER.info("Loaded {} germline mutations", sampleId);

        return variants;
    }
}
