package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS_DELIM;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.parseAltTranscriptInfo;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.driver.panel.ReportablePredicate;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantDrivers
{
    private final Map<VariantType,Integer> mVariantTypeCounts;
    private final Map<VariantType,Integer> mVariantTypeCountsBiallelic;

    private final OncoDrivers mOncoDrivers;
    private final TsgDrivers mTsgDrivers;

    public SomaticVariantDrivers(final DriverGenePanel panel)
    {
        mVariantTypeCounts = Maps.newHashMap();
        mVariantTypeCountsBiallelic = Maps.newHashMap();

        mOncoDrivers = new OncoDrivers(panel);
        mTsgDrivers = new TsgDrivers(panel);
    }

    public void overrideVariantCounts(final Map<VariantType,Integer> variantTypeCounts, final Map<VariantType,Integer> biallelicCounts)
    {
        mVariantTypeCounts.putAll(variantTypeCounts);
        mVariantTypeCountsBiallelic.putAll(biallelicCounts);
    }

    public boolean checkSomaticVariant(final SomaticVariant variant)
    {
        // return true if the variant is a reportable TSG or onocogene
        if(variant.isFiltered())
            return false;

        mVariantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);

        if(variant.biallelic())
        {
            mVariantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);
        }

        if(mOncoDrivers.checkVariant(variant))
            return true;

        if(mTsgDrivers.checkVariant(variant))
            return true;

        return false;
    }

    public void addPhasedReportableVariant(final SomaticVariant variant, final SomaticVariant phasedReportedVariant)
    {
        // add the phased unreported variant to the same driver cache as its reported variant
        if(mOncoDrivers.hasVariant(phasedReportedVariant))
            mOncoDrivers.addVariant(variant);

        if(mTsgDrivers.hasVariant(phasedReportedVariant))
            mTsgDrivers.addVariant(variant);
    }

    protected static boolean isReportable(final ReportablePredicate predicate, final SomaticVariant variant)
    {
        return predicate.isReportable(variant.variantImpact(), variant.type(), variant.isHotspot());
    }

    public static boolean hasTranscriptCodingEffect(final VariantImpact variantImpact, final VariantType variantType, final String transcript)
    {
        if(variantImpact.CanonicalTranscript.equals(transcript))
            return hasCodingEffect(variantType, variantImpact.CanonicalCodingEffect);

        if(!variantImpact.OtherReportableEffects.isEmpty())
        {
            List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variantImpact.OtherReportableEffects);
            return altTransEffects.stream().filter(x -> x.TransName.equals(transcript)).anyMatch(x -> hasCodingEffect(variantType, x.Effect));
        }

        return false;
    }

    public static CodingEffect getWorstReportableCodingEffect(final VariantImpact variantImpact)
    {
        CodingEffect topCodingEffect = variantImpact.CanonicalCodingEffect;

        if(!variantImpact.OtherReportableEffects.isEmpty())
        {
            List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variantImpact.OtherReportableEffects);

            for(AltTranscriptReportableInfo altTransInfo : altTransEffects)
            {
                if(altTransInfo.Effect.ordinal() < topCodingEffect.ordinal())
                {
                    topCodingEffect = altTransInfo.Effect;
                }
            }
        }

        return topCodingEffect;
    }

    public static boolean hasCodingEffect(final VariantType type, final CodingEffect codingEffect)
    {
        DriverImpact impact = DriverImpact.select(type, codingEffect);
        return impact != DriverImpact.UNKNOWN;
    }

    public static void addReportableTranscriptList(
            final VariantType variantType, final VariantContext variantContext, final VariantImpact variantImpact)
    {
        // check alt transcript status vs canonical
        if(variantImpact.OtherReportableEffects.isEmpty())
            return;

        StringJoiner reportableTranscripts = new StringJoiner(REPORTABLE_TRANSCRIPTS_DELIM);

        boolean canonicalReportable = hasCodingEffect(variantType, variantImpact.CanonicalCodingEffect);

        if(canonicalReportable)
            reportableTranscripts.add(variantImpact.CanonicalTranscript);

        List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variantImpact.OtherReportableEffects);

        for(AltTranscriptReportableInfo altTransEffect : altTransEffects)
        {
            if(hasCodingEffect(variantType, altTransEffect.Effect))
                reportableTranscripts.add(altTransEffect.TransName);
        }

        variantContext.getCommonInfo().putAttribute(REPORTABLE_TRANSCRIPTS, reportableTranscripts.toString());
    }

    public List<DriverCatalog> buildCatalog(final Map<String,List<GeneCopyNumber>> geneCopyNumberMap, final List<DriverSourceData> driverSourceData)
    {
        final List<DriverCatalog> result = Lists.newArrayList();

        result.addAll(mOncoDrivers.findDrivers(geneCopyNumberMap, mVariantTypeCounts, mVariantTypeCountsBiallelic, driverSourceData));
        result.addAll(mTsgDrivers.findDrivers(geneCopyNumberMap, mVariantTypeCounts, mVariantTypeCountsBiallelic, driverSourceData));

        return result;
    }

    protected static Map<DriverImpact,Integer> groupByImpact(final List<SomaticVariant> variants)
    {
        Map<DriverImpact,Integer> map = Maps.newHashMap();

        for(SomaticVariant variant : variants)
        {
            DriverImpact driverImpact = DriverImpact.select(variant.type(), variant.variantImpact().CanonicalCodingEffect);
            Integer count = map.get(driverImpact);
            map.put(driverImpact, count != null ? count + 1 : 1);
        }

        return map;
    }
}
