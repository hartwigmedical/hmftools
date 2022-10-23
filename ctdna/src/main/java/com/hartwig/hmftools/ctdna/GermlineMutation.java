package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.CategoryType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.ctdna.PointMutation.generateMutationSequence;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.VariantSelection.addRegisteredLocation;
import static com.hartwig.hmftools.ctdna.VariantSelection.isNearRegisteredLocation;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;

import org.apache.commons.compress.utils.Lists;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class GermlineMutation implements Variant
{
    private final GermlineVariant mVariant;
    private String mSequence;

    public GermlineMutation(final GermlineVariant variant)
    {
        mVariant = variant;

        mSequence = "";
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
    public String sequence() { return mSequence; }

    @Override
    public double copyNumber() { return mVariant.adjustedCopyNumber(); }

    @Override
    public double vaf() { return mVariant.adjustedVAF(); }

    @Override
    public double gc() { return VariantUtils.calcGcPercent(mSequence); }

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
        mSequence = generateMutationSequence(
                refGenome, config, mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
    }

    @Override
    public boolean checkAndRegisterLocation(final Map<String,List<Integer>> registeredLocations)
    {
        if(isNearRegisteredLocation(registeredLocations, mVariant.chromosome(), mVariant.position()))
            return false;

        addRegisteredLocation(registeredLocations, mVariant.chromosome(), mVariant.position());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<Variant> loadGermlineMutations(final String sampleId, final PvConfig config)
    {
        List<Variant> variants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(config.PurpleDir, sampleId);

        try
        {
            List<GermlineVariant> germlineVariants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);

            germlineVariants.stream().filter(x -> x.reported()).forEach(x -> variants.add(new GermlineMutation(x)));

            PV_LOGGER.debug("sample({}) loaded {} germline variants", sampleId, variants.size());
        }
        catch(Exception e)
        {
            PV_LOGGER.error("failed to read germline VCF file({}): {}", vcfFile, e.toString());
        }

        return variants;
    }
}
