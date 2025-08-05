package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_CLONAL_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.SUBCLONAL_MUTATION;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.MAX_INDEL_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.MAX_INSERT_BASES;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.REPEAT_COUNT_MAX_LOWER;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.SUBCLONAL_LIKELIHOOD_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.VAF_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeGenerator.generateMutationProbe;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.Strings;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.wisp.CategoryType;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.ProbeFactory;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SomaticMutation extends Variant
{
    private final VariantContextDecorator mVariantDecorator;
    private final int mTumorDepth;
    private final double mTumorAF;
    private final int mLocationHash;

    private static final Logger LOGGER = LogManager.getLogger(SomaticMutation.class);

    public SomaticMutation(final VariantContext variantContext, final String sampleId)
    {
        mVariantDecorator = new VariantContextDecorator(variantContext);

        Genotype genotype = variantContext.getGenotype(sampleId);

        if(genotype != null)
        {
            mTumorDepth = genotype.getAD()[1];
            mTumorAF = mVariantDecorator.allelicDepth(sampleId).alleleFrequency();
        }
        else
        {
            mTumorDepth = 0;
            mTumorAF = mVariantDecorator.adjustedVaf();
        }

        mLocationHash = Integer.parseInt(Strings.reverseString(String.valueOf(variantContext.getStart())));
    }

    public int locationHash()
    {
        return mLocationHash;
    }

    @Override
    public CategoryType categoryType()
    {
        if(mVariantDecorator.reported())
        {
            return REPORTABLE_MUTATION;
        }

        boolean isCoding = mVariantDecorator.variantImpact().CanonicalCodingEffect != CodingEffect.NONE
                && mVariantDecorator.variantImpact().CanonicalCodingEffect != CodingEffect.UNDEFINED;

        boolean isSubclonal = subclonalLikelihood() >= SUBCLONAL_LIKELIHOOD_MIN;

        if(mVariantDecorator.type() == SNP)
        {
            if(isCoding)
            {
                return OTHER_CODING_MUTATION;
            }

            if(!isSubclonal)
            {
                return OTHER_CLONAL_MUTATION;
            }
        }

        // unused for now
        // if(isSubclonal)
        //    return SUBCLONAL_MUTATION;

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
        return mVariantDecorator.variantImpact().GeneName;
    }

    @Override
    public double copyNumber()
    {
        return mVariantDecorator.adjustedCopyNumber();
    }

    @Override
    public double vaf()
    {
        return mTumorAF;
    }

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
    public int tumorFragments()
    {
        return mTumorDepth;
    }

    @Override
    public boolean hasPhaseVariants()
    {
        return mVariantDecorator.context().hasAttribute(LOCAL_PHASE_SET);
    }

    @Override
    public boolean reported()
    {
        return mVariantDecorator.reported();
    }

    @Override
    public void generateProbe(final RefGenomeInterface refGenome, final ProbeFactory probeFactory, final PanelCoverage coverage)
    {
        ProbeGenerationResult result = generateMutationProbe(
                mVariantDecorator.chromosome(), mVariantDecorator.position(), mVariantDecorator.ref(), mVariantDecorator.alt(),
                targetMetadata(), refGenome, probeFactory, coverage);
        setProbeGenResult(result);
    }

    private double subclonalLikelihood()
    {
        return mVariantDecorator.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);
    }

    @Override
    public boolean checkFilters()
    {
        return !reported();
    }

    @Override
    public boolean passNonReportableFilters(boolean useLowerLimits)
    {
        if(!passesGcRatioLimit(probe().gcContent(), useLowerLimits))
        {
            return false;
        }

        if(categoryType() != SUBCLONAL_MUTATION && vaf() < VAF_MIN)
        {
            return false;
        }

        if(!passesFragmentCountLimit(tumorFragments(), useLowerLimits))
        {
            return false;
        }

        int repeatCountMax = max(
                mVariantDecorator.repeatCount(), mVariantDecorator.context().getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0));

        int maxRepeatCount = useLowerLimits ? REPEAT_COUNT_MAX_LOWER : REPEAT_COUNT_MAX;
        if(repeatCountMax > maxRepeatCount)
        {
            return false;
        }

        if(mVariantDecorator.type() == VariantType.INDEL)
        {
            if(max(mVariantDecorator.alt().length(), mVariantDecorator.ref().length()) > MAX_INDEL_LENGTH)
            {
                return false;
            }
        }

        GermlineStatus germlineStatus = GermlineStatus.valueOf(
                mVariantDecorator.context().getAttributeAsString(PURPLE_GERMLINE_INFO, GermlineStatus.UNKNOWN.toString()));

        if(germlineStatus == GermlineStatus.AMPLIFICATION || germlineStatus == GermlineStatus.NOISE)
        {
            return false;
        }

        return true;
    }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariantDecorator.chromosome(), mVariantDecorator.position()))
        {
            return false;
        }

        registeredLocations.addRegisteredLocation(mVariantDecorator.chromosome(), mVariantDecorator.position());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<Variant> loadSomatics(final String sampleId, final String purpleDir)
    {
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);

        List<Variant> variants = Lists.newArrayList();

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            String error = "failed to read somatic vcf: " + vcfFile;
            LOGGER.error(error);
            throw new RuntimeException(error);
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
            {
                continue;
            }

            String alt = VariantContextDecorator.getAlt(variantContext);
            if(alt.length() >= MAX_INSERT_BASES)
            {
                continue;
            }

            variants.add(new SomaticMutation(variantContext, sampleId));
        }

        LOGGER.info("loaded {} somatic variants from vcf({})", variants.size(), vcfFile);

        return variants;
    }
}
