package com.hartwig.hmftools.common.variant.impact;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class AltTranscriptReportableInfo
{
    public final String TransName;
    public final String HgvsCoding;
    public final String HgvsProtein;
    public final String Effects;
    public final CodingEffect Effect;

    public static final String VAR_IMPACT_OTHER_REPORT_ITEM_DELIM = "|";
    public static final String VAR_IMPACT_OTHER_REPORT_DELIM = "-";
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
}
