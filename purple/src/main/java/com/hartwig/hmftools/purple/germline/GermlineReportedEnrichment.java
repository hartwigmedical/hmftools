package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.addReportableTranscriptList;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment
{
    private final Map<String,DriverGene> mDriverGeneMap;
    private final List<GermlineVariant> mBuffer;

    private final Set<String> mSomaticReportedGenes;
    private final List<GermlineVariant> mCandidateReportableVariants;

    private static final double MIN_VARIANT_COPY_NUMBER = 0.5;

    public GermlineReportedEnrichment(final Map<String,DriverGene> driverGenes, final Set<String> somaticReportedGenes)
    {
        mDriverGeneMap = driverGenes;
        mSomaticReportedGenes = somaticReportedGenes;

        mCandidateReportableVariants = Lists.newArrayList();

        mBuffer = Lists.newArrayList();
    }

    public List<GermlineVariant> candidateReportableVariants() { return mCandidateReportableVariants; }

    public void processVariant(final GermlineVariant variant)
    {
        mBuffer.add(variant);
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
    }

    public void flush()
    {
        Map<String,List<GermlineVariant>> geneVariantsMap = Maps.newHashMap();

        for(GermlineVariant variant : mBuffer)
        {
            if(!mDriverGeneMap.containsKey(variant.gene()))
                continue;

            if(isCandidateReportable(variant.decorator()))
            {
                mCandidateReportableVariants.add(variant);

                List<GermlineVariant> geneVariants = geneVariantsMap.get(variant.gene());

                if(geneVariants == null)
                {
                    geneVariants = Lists.newArrayList();
                    geneVariantsMap.put(variant.gene(), geneVariants);
                }

                geneVariants.add(variant);
            }
        }

        // check reportability status
        for(Map.Entry<String,List<GermlineVariant>> entry : geneVariantsMap.entrySet())
        {
            String gene = entry.getKey();
            List<GermlineVariant> geneVariants = entry.getValue();

            DriverGene driverGene = mDriverGeneMap.get(gene);

            boolean hasMultipleUnphasedHits = false;

            for(GermlineVariant variant : geneVariants)
            {
                if(geneVariants.size() > 1)
                {
                    hasMultipleUnphasedHits |= hasUnphasedSameGeneVariant(variant, geneVariants);
                }
            }

            if(driverGene.reportGermlineVariant() != NONE || driverGene.reportGermlineHotspot() != NONE)
            {
                Set<String> multiHitGenes = Sets.newHashSet(mSomaticReportedGenes);

                if(hasMultipleUnphasedHits)
                    multiHitGenes.add(gene);

                for(GermlineVariant variant : geneVariants)
                {
                    if(isReportable(variant, driverGene.reportGermlineHotspot(), driverGene.reportGermlineVariant(), multiHitGenes))
                    {
                        variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
                        addReportableTranscriptList(variant.type(), variant.context(), variant.variantImpact());
                    }
                }
            }
        }

        mBuffer.clear();
    }

    private boolean isCandidateReportable(final VariantContextDecorator variant)
    {
        // check reportable criteria irrespective of the reportable flags set in the driver gene catalog
        if(variant.gene().isEmpty())
            return false;

        if(!mDriverGeneMap.containsKey(variant.gene()))
            return false;

        if(!variant.isPass())
            return false;

        return variant.isGermlinePathogenic() || variant.variantImpact().CanonicalSpliceRegion;
    }

    private static boolean isReportable(
            final GermlineVariant variant, final DriverGeneGermlineReporting hotspotReporting,
            final DriverGeneGermlineReporting variantReporting, final Set<String> genesWithMultipleUnphasedHits)
    {
        // check the driver gene configuration
        if(!variant.decorator().isGermlinePathogenic())
            return false;

        DriverGeneGermlineReporting reporting = variant.isHotspot() ? hotspotReporting : variantReporting;

        if(reporting == NONE)
            return false;

        if(reporting == ANY)
            return true;

        if(isVariantLost(variant, MIN_VARIANT_COPY_NUMBER))
            return false;

        if(reporting == VARIANT_NOT_LOST)
            return true;

        if(reporting == WILDTYPE_LOST)
            return variant.biallelic() || genesWithMultipleUnphasedHits.contains(variant.gene());

        return false;
    }

    private static boolean isVariantLost(final GermlineVariant variant, double minVariantCopyNumber)
    {
        return Doubles.lessThan(variant.decorator().variantCopyNumber(), minVariantCopyNumber);
    }

    private static boolean hasUnphasedSameGeneVariant(final GermlineVariant variant, final List<GermlineVariant> otherVariants)
    {
        for(GermlineVariant otherVariant : otherVariants)
        {
            if(otherVariant == variant)
                continue;

            if(!otherVariant.gene().equals(variant.gene()))
                continue;

            if(variant.decorator().localPhaseSet() == null
            || otherVariant.decorator().localPhaseSet() == null
            || !Objects.equals(variant.decorator().localPhaseSet(), otherVariant.decorator().localPhaseSet()))
            {
                return true;
            }
        }

        return false;
    }
}
