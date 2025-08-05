package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeGenerator.generateMutationProbe;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.common.wisp.CategoryType;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.ProbeFactory;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;

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
    public void generateProbe(final RefGenomeInterface refGenome, final ProbeFactory probeFactory, final PanelCoverage coverage)
    {
        ProbeGenerationResult result = generateMutationProbe(
                mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(),
                targetMetadata(), refGenome, probeFactory, coverage);
        setProbeGenResult(result);
    }

    @Override
    public boolean checkFilters()
    {
        return false;
    }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
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

    public static List<Variant> loadGermlineMutations(final String sampleId, final String purpleDir)
    {
        List<Variant> variants = Lists.newArrayList();

        // TODO: unused - bug?
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);

        // TODO: better error handling?
        try
        {
            List<GermlineVariant> germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);

            germlineVariants.stream().filter(GermlineVariant::reported).forEach(x -> variants.add(new GermlineMutation(x)));

            LOGGER.debug("sample({}) loaded {} germline variants", sampleId, variants.size());
        }
        catch(Exception e)
        {
            LOGGER.error("sample({}) failed to load germline variants from file: {}", sampleId, vcfFile);
        }

        return variants;
    }
}
