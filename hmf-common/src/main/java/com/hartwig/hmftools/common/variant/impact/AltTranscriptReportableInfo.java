package com.hartwig.hmftools.common.variant.impact;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.util.Strings;

public class AltTranscriptReportableInfo
{
    public final String GeneName;
    public final String TransName;
    public final String HgvsCoding;
    public final String HgvsProtein;
    public final String Effects;
    public final CodingEffect Effect;

    public static final String VAR_IMPACT_OTHER_REPORT_DELIM = "--"; // a single hyphen conflicts with the HGVS coding annotation
    private static final String VAR_IMPACT_OTHER_REPORT_ITEM_DELIM = "|";
    private static final int VAR_IMPACT_OTHER_REPORT_ITEM_COUNT = 6;
    private static final int VAR_IMPACT_OTHER_REPORT_ITEM_COUNT_OLD = 5;

    public AltTranscriptReportableInfo(
            final String geneName, final String transName, final String hgvsCoding, final String hgvsProtein, final String effects, final CodingEffect codingEffect)
    {
        GeneName = geneName;
        TransName = transName;
        HgvsCoding = hgvsCoding;
        HgvsProtein = hgvsProtein;
        Effects = effects;
        Effect = codingEffect;
    }

    public static List<AltTranscriptReportableInfo> parseAltTranscriptInfo(final String otherReportableEffects)
    {
        List<AltTranscriptReportableInfo> altTransInfos = new ArrayList<>();

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
        if(transValues.length < VAR_IMPACT_OTHER_REPORT_ITEM_COUNT_OLD)
            return null;

        int itemIndex = 0;
        String geneName = "";

        if(transValues.length == VAR_IMPACT_OTHER_REPORT_ITEM_COUNT)
        {
            geneName = transValues[itemIndex++];
        }

        return new AltTranscriptReportableInfo(
                geneName,
                transValues[itemIndex++],
                transValues[itemIndex++],
                transValues[itemIndex++],
                transValues[itemIndex++],
                CodingEffect.valueOf(transValues[itemIndex]));
    }

    public String serialise() { return serialise(GeneName, TransName, HgvsCoding, HgvsProtein, Effects, Effect); }

    public static String serialise(
            final String geneName, final String transName, final String hgvsCoding, final String hgvsProtein, final String effects, final CodingEffect codingEffect)
    {
        // eg CDKN2A|ENST00000579755|c.209_210delCCinsTT|p.Pro70Leu|missense_variant|MISSENSE;
        StringJoiner sj = new StringJoiner(VAR_IMPACT_OTHER_REPORT_ITEM_DELIM);
        sj.add(geneName);
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

