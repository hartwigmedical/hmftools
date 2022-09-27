package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.ctdna.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.ctdna.PvConfig.MAX_INSERT_BASES;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantVcfTags;

import org.apache.commons.compress.utils.Lists;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

public class PointMutation implements Variant
{
    private final VariantContextDecorator mVariantDecorator;
    private final int mTumorDepth;
    private String mSequence;

    public PointMutation(final VariantContext variantContext, final String sampleId)
    {
        mVariantDecorator = new VariantContextDecorator(variantContext);

        final Genotype genotype = variantContext.getGenotype(sampleId);

        if(genotype != null)
            mTumorDepth = genotype.getAD()[1];
        else
            mTumorDepth = 0;
        mSequence = "";
    }

    @Override
    public CategoryType categoryType()
    {
        if(mVariantDecorator.reported())
            return REPORTABLE_MUTATION;

        if(mVariantDecorator.variantImpact().CanonicalCodingEffect != CodingEffect.NONE)
            return OTHER_CODING_MUTATION;

        return OTHER_MUTATION;
    }

    @Override
    public String description()
    {
        return format("%s:%s %s>%s %s", mVariantDecorator.chromosome(), mVariantDecorator.position(),
                mVariantDecorator.ref(), mVariantDecorator.alt(), mVariantDecorator.type());
    }

    @Override
    public String gene()
    {
        return mVariantDecorator.variantImpact().CanonicalGeneName;
    }

    @Override
    public String sequence() { return mSequence; }

    @Override
    public double copyNumber() { return mVariantDecorator.adjustedCopyNumber(); }

    @Override
    public double vaf() { return mVariantDecorator.adjustedVaf(); }

    @Override
    public int tumorFragments() { return mTumorDepth; }

    @Override
    public boolean hasPhaseVariants()
    {
        return mVariantDecorator.context().hasAttribute(VariantVcfTags.LOCAL_PHASE_SET);
    }

    @Override
    public boolean reported() { return mVariantDecorator.reported(); }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final PvConfig config)
    {
        int variantPosition = mVariantDecorator.position();
        String alt = mVariantDecorator.alt();
        int altLength = alt.length();
        int refLength = mVariantDecorator.ref().length();
        int startLength = config.ProbeLength / 2 - altLength / 2;
        int startPos = mVariantDecorator.position() - startLength;

        String basesStart = refGenome.getBaseString(mVariantDecorator.chromosome(), startPos, variantPosition - 1);
        int endBaseLength = config.ProbeLength - basesStart.length() - altLength;

        int postPosition = variantPosition + refLength;
        String basesEnd = refGenome.getBaseString(mVariantDecorator.chromosome(), postPosition, postPosition + endBaseLength - 1);

        mSequence = basesStart + alt + basesEnd;

        if(mSequence.length() != config.ProbeLength)
        {
            PV_LOGGER.error("variant({}) invalid sequenceLength({}): {}", description(), mSequence.length(), mSequence);
        }
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<Variant> loadSomatics(final PvConfig config)
    {
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(config.PurpleDir, config.SampleId);

        List<Variant> variants = Lists.newArrayList();

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            CompoundFilter filter = new CompoundFilter(true);
            filter.add(new PassingVariantFilter());

            for(VariantContext variantContext : reader.iterator())
            {
                if(!filter.test(variantContext))
                    continue;

                String alt = VariantContextDecorator.getAlt(variantContext);
                if(alt.length() >= MAX_INSERT_BASES)
                    continue;

                variants.add(new PointMutation(variantContext, config.SampleId));
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read somatic VCF file({}): {}", vcfFile, e.toString());
            return variants;
        }

        PV_LOGGER.info("loaded {} somatic variants from vcf({})", variants.size(), vcfFile);

        return variants;
    }
}
