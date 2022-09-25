package com.hartwig.hmftools.common.variant.impact;

import java.util.List;
import java.util.StringJoiner;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class AltTranscriptReportableInfo
{
    public final String TransName;
    public final String HgvsCoding;
    public final String HgvsProtein;
    public final String Effects;
    public final CodingEffect Effect;

    public static final String VAR_IMPACT_OTHER_REPORT_ITEM_DELIM = "|";
    public static final String VAR_IMPACT_OTHER_REPORT_DELIM = "--"; // a single hyphen conflicts with the HGVS coding annotation
    public static final int VAR_IMPACT_OTHER_REPORT_ITEM_COUNT = 5;

    public AltTranscriptReportableInfo(
            final String transName, final String hgvsCoding, final String hgvsProtein, final String effects, final CodingEffect codingEffect)
    {
        TransName = transName;
        HgvsCoding = hgvsCoding;
        HgvsProtein = hgvsProtein;
        Effects = effects;
        Effect = codingEffect;
    }

    public static List<AltTranscriptReportableInfo> parseAltTranscriptInfo(final String otherReportableEffects)
    {
        List<AltTranscriptReportableInfo> altTransInfos = Lists.newArrayList();

        String[] altTransInfoItems = otherReportableEffects.split(VAR_IMPACT_OTHER_REPORT_DELIM, -1);

        for(String altTransInfoStr : altTransInfoItems)
        {
            AltTranscriptReportableInfo altTransInfo = parse(altTransInfoStr);

            if(altTransInfo != null)
                altTransInfos.add(altTransInfo);
        }

        return altTransInfos;
    }

    public static AltTranscriptReportableInfo parse(final String transInfo)
    {
        String[] transValues = transInfo.split("\\" + VAR_IMPACT_OTHER_REPORT_ITEM_DELIM, -1);
        if(transValues.length != VAR_IMPACT_OTHER_REPORT_ITEM_COUNT)
            return null;

        return new AltTranscriptReportableInfo(
                transValues[0], transValues[1], transValues[2], transValues[3], CodingEffect.valueOf(transValues[4]));
    }

    public String serialise() { return serialise(TransName, HgvsCoding, HgvsProtein, Effects, Effect); }

    public static String serialise(
            final String transName, final String hgvsCoding, final String hgvsProtein, final String effects, final CodingEffect codingEffect)
    {
        // eg ENST00000579755|c.209_210delCCinsTT|p.Pro70Leu|missense_variant|MISSENSE;
        StringJoiner sj = new StringJoiner(VAR_IMPACT_OTHER_REPORT_ITEM_DELIM);
        sj.add(transName);
        sj.add(hgvsCoding);
        sj.add(hgvsProtein);
        sj.add(effects);
        sj.add(codingEffect.toString());
        return sj.toString();
    }

    // convenience methods for Protect

    @Nullable
    public static AltTranscriptReportableInfo getFirstAltTranscriptInfo(final String otherReportableEffects)
    {
        if(otherReportableEffects.isEmpty())
            return null;

        return parseAltTranscriptInfo(otherReportableEffects).get(0);
    }

    public static String firstOtherTranscript(final String otherReportedEffects)
    {
        if(otherReportedEffects.isEmpty())
            return Strings.EMPTY;

        return getFirstAltTranscriptInfo(otherReportedEffects).TransName;
    }

    public static String firstOtherEffects(final String otherReportedEffects)
    {
        if(otherReportedEffects.isEmpty())
            return Strings.EMPTY;

        return getFirstAltTranscriptInfo(otherReportedEffects).Effects;
    }

    public static String firstOtherHgvsCodingImpact(final String otherReportedEffects)
    {
        if(otherReportedEffects.isEmpty())
            return Strings.EMPTY;

        return getFirstAltTranscriptInfo(otherReportedEffects).HgvsCoding;
    }

    public static String firstOtherHgvsProteinImpact(final String otherReportedEffects)
    {
        if(otherReportedEffects.isEmpty())
            return Strings.EMPTY;

        return getFirstAltTranscriptInfo(otherReportedEffects).HgvsProtein;
    }

    public static CodingEffect firstOtherCodingEffect(final String otherReportedEffects)
    {
        if(otherReportedEffects.isEmpty())
            return CodingEffect.UNDEFINED;

        return getFirstAltTranscriptInfo(otherReportedEffects).Effect;
    }

}

