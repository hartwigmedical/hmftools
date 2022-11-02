package com.hartwig.hmftools.ctdna;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.ctdna.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.ctdna.CategoryType.SUBCLONAL_MUTATION;
import static com.hartwig.hmftools.ctdna.PvConfig.DEFAULT_GC_THRESHOLD_MAX;
import static com.hartwig.hmftools.ctdna.PvConfig.DEFAULT_GC_THRESHOLD_MIN;
import static com.hartwig.hmftools.ctdna.PvConfig.DEFAULT_MAPPABILITY_MIN;
import static com.hartwig.hmftools.ctdna.PvConfig.DEFAULT_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.ctdna.PvConfig.DEFAULT_SUBCLONAL_LIKELIHOOD_MIN;
import static com.hartwig.hmftools.ctdna.PvConfig.MAX_INDEL_LENGTH;
import static com.hartwig.hmftools.ctdna.PvConfig.MAX_INSERT_BASES;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.VariantSelection.addRegisteredLocation;
import static com.hartwig.hmftools.ctdna.VariantSelection.isNearRegisteredLocation;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.Strings;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.CommonVcfTags;

import org.apache.commons.compress.utils.Lists;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

public class PointMutation extends Variant
{
    private final VariantContextDecorator mVariantDecorator;
    private final int mTumorDepth;
    private String mSequence;
    private int mLocationHash;

    public PointMutation(final VariantContext variantContext, final String sampleId)
    {
        mVariantDecorator = new VariantContextDecorator(variantContext);

        final Genotype genotype = variantContext.getGenotype(sampleId);

        if(genotype != null)
            mTumorDepth = genotype.getAD()[1];
        else
            mTumorDepth = 0;

        mLocationHash = Integer.parseInt(Strings.reverseString(String.valueOf(variantContext.getStart())));

        mSequence = "";
    }

    public VariantContextDecorator variantDecorator() { return mVariantDecorator; }
    public int locationHash() { return mLocationHash; }

    @Override
    public CategoryType categoryType()
    {
        if(mVariantDecorator.reported())
            return REPORTABLE_MUTATION;

        if(subclonalLikelihood() >= DEFAULT_SUBCLONAL_LIKELIHOOD_MIN)
            return SUBCLONAL_MUTATION;

        if(mVariantDecorator.variantImpact().CanonicalCodingEffect != CodingEffect.NONE
        && mVariantDecorator.variantImpact().CanonicalCodingEffect != CodingEffect.UNDEFINED)
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
    public double gc() { return VariantUtils.calcGcPercent(mSequence); }

    @Override
    public String otherData()
    {
        return format("Map=%.2f Repeats=%d/%d SubClonal=%.2f GermlineStatus=%s",
            mVariantDecorator.mappability(), mVariantDecorator.repeatCount(),
                mVariantDecorator.context().getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0),
                subclonalLikelihood(),
                mVariantDecorator.context().getAttributeAsString(PURPLE_GERMLINE_INFO, ""));
    }

    @Override
    public int tumorFragments() { return mTumorDepth; }

    @Override
    public boolean hasPhaseVariants()
    {
        return mVariantDecorator.context().hasAttribute(LOCAL_PHASE_SET);
    }

    @Override
    public boolean reported() { return mVariantDecorator.reported(); }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final PvConfig config)
    {
        mSequence = generateMutationSequence(
                refGenome, config, mVariantDecorator.chromosome(), mVariantDecorator.position(), mVariantDecorator.ref(),
                mVariantDecorator.alt());
    }

    protected static String generateMutationSequence(
            final RefGenomeInterface refGenome, final PvConfig config,
            final String chromosome, final int position, final String ref, final String alt)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = config.ProbeLength / 2 - altLength / 2;
        int startPos = position - startLength;

        String basesStart = refGenome.getBaseString(chromosome, startPos, position - 1);
        int endBaseLength = config.ProbeLength - basesStart.length() - altLength;

        int postPosition = position + refLength;
        String basesEnd = refGenome.getBaseString(chromosome, postPosition, postPosition + endBaseLength - 1);

        String sequence = basesStart + alt + basesEnd;

        if(sequence.length() != config.ProbeLength)
        {
            PV_LOGGER.error("variant({}:{} {}->{}) invalid sequenceLength({}): {}",
                    chromosome, position, ref, alt, sequence.length(), sequence);
        }

        return sequence;
    }

    private double subclonalLikelihood()
    {
        return mVariantDecorator.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);
    }

    @Override
    public boolean passNonReportableFilters(final PvConfig config)
    {
        if(gc() < DEFAULT_GC_THRESHOLD_MIN || gc() > DEFAULT_GC_THRESHOLD_MAX)
            return false;

        if(categoryType() != SUBCLONAL_MUTATION && vaf() < config.VafMin)
            return false;

        if(tumorFragments() < config.FragmentCountMin)
            return false;

        if(mVariantDecorator.mappability() < DEFAULT_MAPPABILITY_MIN)
            return false;

        int repeatCountMax = max(
                mVariantDecorator.repeatCount(), mVariantDecorator.context().getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0));

        if(repeatCountMax > DEFAULT_REPEAT_COUNT_MAX)
            return false;

        if(mVariantDecorator.type() == VariantType.INDEL)
        {
            if(max(mVariantDecorator.alt().length(), mVariantDecorator.ref().length()) > MAX_INDEL_LENGTH)
                return false;
        }

        GermlineStatus germlineStatus = GermlineStatus.valueOf(
                mVariantDecorator.context().getAttributeAsString(PURPLE_GERMLINE_INFO, GermlineStatus.UNKNOWN.toString()));

        if(germlineStatus == GermlineStatus.AMPLIFICATION || germlineStatus == GermlineStatus.NOISE)
            return false;

        return true;
    }


    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariantDecorator.chromosome(), mVariantDecorator.position()))
            return false;

        registeredLocations.addRegisteredLocation(mVariantDecorator.chromosome(), mVariantDecorator.position());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<Variant> loadSomatics(final String sampleId, final PvConfig config) throws Exception
    {
        String purpleDir = PvConfig.getSampleFilePath(sampleId, config.PurpleDir);
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);

        List<Variant> variants = Lists.newArrayList();

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

            variants.add(new PointMutation(variantContext, sampleId));
        }

        PV_LOGGER.info("loaded {} somatic variants from vcf({})", variants.size(), vcfFile);

        return variants;
    }
}
