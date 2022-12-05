package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.CategoryType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.ctdna.PointMutation.generateMutationSequence;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.VariantSelection.addRegisteredLocation;
import static com.hartwig.hmftools.ctdna.VariantSelection.isNearRegisteredLocation;

import java.util.List;
import java.util.Map;

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
    public CategoryType categoryType() { return GERMLINE_MUTATION; }

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
    public int tumorFragments() { return mVariant.alleleReadCount(); }

    @Override
    public boolean hasPhaseVariants()
    {
        return false;
    }

    @Override
    public boolean reported() { return true; }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final PvConfig config)
    {
        String sequence = generateMutationSequence(
                refGenome, config, mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
        setSequence(sequence);
    }

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

    public static List<Variant> loadGermlineMutations(final String sampleId, final PvConfig config) throws Exception
    {
        List<Variant> variants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String purpleDir = PvConfig.getSampleFilePath(sampleId, config.PurpleDir);
        String vcfFile = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);

        List<GermlineVariant> germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);

        germlineVariants.stream().filter(x -> x.reported()).forEach(x -> variants.add(new GermlineMutation(x)));

        PV_LOGGER.debug("sample({}) loaded {} germline variants", sampleId, variants.size());

        return variants;
    }
}
