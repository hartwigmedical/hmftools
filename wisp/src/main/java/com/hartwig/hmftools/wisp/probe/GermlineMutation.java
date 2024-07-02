package com.hartwig.hmftools.wisp.probe;

import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.generateMutationSequence;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class GermlineMutation extends Variant
{
    private final GermlineVariant mVariant;

    public GermlineMutation(final GermlineVariant variant)
    {
        mVariant = variant;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.GERMLINE_MUTATION; }

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
    public double copyNumber() { return mVariant.adjustedCopyNumber(); }

    @Override
    public double vaf() { return mVariant.adjustedVAF(); }

    @Override
    public int tumorFragments() { return mVariant.allelicDepth().AlleleReadCount; }

    @Override
    public boolean hasPhaseVariants()
    {
        return false;
    }

    @Override
    public boolean reported() { return true; }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final ProbeConfig config)
    {
        String sequence = generateMutationSequence(
                refGenome, config.ProbeLength, mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
        setSequence(sequence);
    }

    @Override
    boolean checkFilters() { return false; }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariant.chromosome(), mVariant.position()))
            return false;

        registeredLocations.addRegisteredLocation(mVariant.chromosome(), mVariant.position());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<Variant> loadGermlineMutations(final String sampleId, final ProbeConfig config)
    {
        List<Variant> variants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String purpleDir = ProbeConfig.getSampleFilePath(sampleId, config.PurpleDir);
        String vcfFile = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);

        try
        {
            List<GermlineVariant> germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);

            germlineVariants.stream().filter(x -> x.reported()).forEach(x -> variants.add(new GermlineMutation(x)));

            CT_LOGGER.debug("sample({}) loaded {} germline variants", sampleId, variants.size());
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load germline variants from file: {}", sampleId, vcfFile);
        }

        return variants;
    }
}
