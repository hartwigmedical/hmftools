package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.parseAltTranscriptInfo;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_ACCEPTOR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.effectsToList;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

public class ReportablePredicate
{
    private final Map<String, DriverGene> mDriverGeneMap;

    public ReportablePredicate(final DriverCategory type, final List<DriverGene> driverGenes)
    {
        mDriverGeneMap = driverGenes.stream()
                .filter(x -> x.likelihoodType().equals(type) && x.reportSomatic())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    public boolean isReportable(final VariantImpact variantImpact, final VariantType type, boolean isHotspot)
    {
        if(isReportable(
            variantImpact.GeneName, type, isHotspot, variantImpact.CanonicalCodingEffect, variantImpact.CanonicalEffect))
        {
            return true;
        }

        if(variantImpact.OtherReportableEffects.isEmpty())
            return false;

        List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variantImpact.OtherReportableEffects);

        return altTransEffects.stream().anyMatch(x ->
                isReportable(variantImpact.GeneName, type, isHotspot, x.Effect, x.Effects));
    }

    public boolean isReportable(
            final String gene, final VariantType type, boolean isHotspot,
            final CodingEffect codingEffect, final String effectsStr)
    {
        DriverGene driverGene = mDriverGeneMap.get(gene);

        if(driverGene == null)
            return false;

        if(isHotspot && driverGene.reportSomaticHotspot())
            return true;

        // rather than take the top driver impact and check its reportability only, check if any present impact is reportable
        // to avoid not reporting a variant with a lower reportable effect with a higher non-reportable one

        List<VariantEffect> variantEffects = effectsToList(effectsStr);

        if(driverGene.reportMissenseAndInframe())
        {
            if(variantEffects.stream().anyMatch(x -> CodingEffect.effect(x) == CodingEffect.MISSENSE))
                return true;
        }

        if(driverGene.reportSplice())
        {
            if(variantEffects.stream().anyMatch(x -> CodingEffect.effect(x) == CodingEffect.SPLICE))
                return true;
        }

        if(driverGene.reportNonsenseAndFrameshift())
        {
            if(variantEffects.stream().anyMatch(x -> CodingEffect.effect(x) == CodingEffect.NONSENSE_OR_FRAMESHIFT))
                return true;
        }

        return false;
    }
}
