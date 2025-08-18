package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SUBCLONAL_LIKELIHOOD_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.util.List;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.panelbuilder.ProbeTarget;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

// Somatic SNV or INDEL.
public class SomaticMutation implements Variant
{
    private final VariantContextDecorator mVariantDecorator;
    private final int mTumorDepth;
    private final double mTumorAF;

    private static final Logger LOGGER = LogManager.getLogger(SomaticMutation.class);

    private SomaticMutation(final VariantContext variantContext, final String sampleId)
    {
        mVariantDecorator = new VariantContextDecorator(variantContext);

        Genotype genotype = variantContext.getGenotype(sampleId);
        if(genotype == null)
        {
            mTumorDepth = 0;
            mTumorAF = mVariantDecorator.adjustedVaf();
        }
        else
        {
            mTumorDepth = genotype.getAD()[1];
            mTumorAF = mVariantDecorator.allelicDepth(sampleId).alleleFrequency();
        }
    }

    public double vaf()
    {
        return mTumorAF;
    }

    public int tumorFragments()
    {
        return mTumorDepth;
    }

    public int repeatCount()
    {
        return max(mVariantDecorator.repeatCount(), mVariantDecorator.context().getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0));
    }

    public int indelLength()
    {
        return abs(mVariantDecorator.alt().length() - mVariantDecorator.ref().length());
    }

    public GermlineStatus germlineStatus()
    {
        return GermlineStatus.valueOf(mVariantDecorator.context()
                .getAttributeAsString(PURPLE_GERMLINE_INFO, GermlineStatus.UNKNOWN.toString()));
    }

    public boolean isCoding()
    {
        return switch(mVariantDecorator.canonicalCodingEffect())
        {
            case NONE, UNDEFINED -> false;
            default -> true;
        };
    }

    public boolean isClonal()
    {
        double subclonalLikelihood = mVariantDecorator.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);
        boolean isSubclonal = subclonalLikelihood >= SAMPLE_SUBCLONAL_LIKELIHOOD_MIN;
        return !isSubclonal;
    }

    // TODO: use in extra info
    private String gene()
    {
        return mVariantDecorator.variantImpact().GeneName;
    }

    @Override
    public boolean isDriver()
    {
        return mVariantDecorator.reported();
    }

    @Override
    public ProbeTarget generateProbeTarget()
    {
        return buildMutationProbe(
                mVariantDecorator.chromosome(), mVariantDecorator.position(), mVariantDecorator.ref(), mVariantDecorator.alt(),
                PROBE_LENGTH);
    }

    @Override
    public TargetMetadata.Type targetType()
    {
        if(isDriver())
        {
            return TargetMetadata.Type.SAMPLE_SNV_INDEL_DRIVER;
        }
        else
        {
            return TargetMetadata.Type.SAMPLE_SNV_INDEL_OTHER;
        }
    }

    @Override
    public String toString()
    {
        return format("%s:%s %s>%s %s",
                mVariantDecorator.chromosome(), mVariantDecorator.position(), mVariantDecorator.ref(), mVariantDecorator.alt(),
                mVariantDecorator.type());
    }

    // Hash code for "randomly" ordering variants while maintaining determinism for comparison purposes.
    public int deterministicHash()
    {
        return mVariantDecorator.chromosome().hashCode() ^ mVariantDecorator.position();
    }

    public static List<SomaticMutation> load(final String sampleId, final String purpleDir)
    {
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);
        if(!vcfFileReader.fileValid())
        {
            throw new RuntimeException("Failed to read somatic vcf: " + vcfFile);
        }

        List<SomaticMutation> variants;
        try(CloseableTribbleIterator<VariantContext> iterator = vcfFileReader.iterator())
        {
            variants = iterator.stream()
                    .filter(variant -> !variant.isFiltered())
                    .map(variant -> new SomaticMutation(variant, sampleId))
                    .toList();
        }

        LOGGER.info("Loaded {} somatic mutations", variants.size());
        variants.forEach(variant -> LOGGER.trace("SomaticMutation: {}", variant));

        return variants;
    }
}
