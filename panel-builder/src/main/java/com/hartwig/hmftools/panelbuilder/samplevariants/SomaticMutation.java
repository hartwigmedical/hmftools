package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_MAX_INDEL;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

// Somatic SNV or INDEL.
public class SomaticMutation extends Variant
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

    @Override
    public boolean isDriver()
    {
        return mVariantDecorator.reported();
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        return buildMutationProbe(
                mVariantDecorator.chromosome(), mVariantDecorator.position(), mVariantDecorator.ref(), mVariantDecorator.alt(),
                PROBE_LENGTH, refGenome);
    }

    // TODO
    public boolean passNonReportableFilters()
    {
        if(!(vaf() >= SAMPLE_VAF_MIN))
        {
            return false;
        }

        if(!(tumorFragments() >= SAMPLE_FRAGMENT_COUNT_MIN))
        {
            return false;
        }

        int repeatCountMax =
                max(mVariantDecorator.repeatCount(), mVariantDecorator.context().getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0));
        if(!(repeatCountMax <= SAMPLE_REPEAT_COUNT_MAX))
        {
            return false;
        }

        if(mVariantDecorator.type() == VariantType.INDEL)
        {
            if(max(mVariantDecorator.alt().length(), mVariantDecorator.ref().length()) > SAMPLE_MAX_INDEL)
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
    public List<ProximateLocations.Location> checkedLocations()
    {
        return List.of(new ProximateLocations.Location(mVariantDecorator.chromosome(), mVariantDecorator.position()));
    }

    @Override
    public String toString()
    {
        return format("%s %s:%s %s>%s",
                mVariantDecorator.type(), mVariantDecorator.chromosome(), mVariantDecorator.position(),
                mVariantDecorator.ref(), mVariantDecorator.alt());
    }

    public static List<SomaticMutation> load(final String sampleId, final String purpleDir)
    {
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);
        if(!vcfFileReader.fileValid())
        {
            throw new RuntimeException("Failed to read somatic vcf: " + vcfFile);
        }

        List<SomaticMutation> variants = vcfFileReader.iterator().stream()
                // TODO: what is this doing?
                .filter(variant -> !variant.isFiltered())
                .map(variant -> new SomaticMutation(variant, sampleId))
                .toList();

        LOGGER.info("Loaded {} somatic mutations", variants.size());
        variants.forEach(variant -> LOGGER.trace("SomaticMutation: {}", variant));

        return variants;
    }
}
