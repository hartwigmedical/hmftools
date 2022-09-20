package com.hartwig.hmftools.purple.drivers;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.commons.compress.utils.Lists;

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

    public static List<AltTranscriptReportableInfo> parseAltTranscriptInfo(final String otherReportableEffects)
    {
        List<AltTranscriptReportableInfo> altTransEffects = Lists.newArrayList();

        String[] otherReportableTranscripts = otherReportableEffects.split(VAR_IMPACT_OTHER_REPORT_DELIM, -1);

        for(String transInfo : otherReportableTranscripts)
        {
            AltTranscriptReportableInfo altTransInfo = AltTranscriptReportableInfo.parse(transInfo);
            if(altTransInfo != null)
                altTransEffects.add(altTransInfo);
        }

        return altTransEffects;
    }

    private static AltTranscriptReportableInfo parse(final String transInfo)
    {
        String[] transValues = transInfo.split("\\" + VAR_IMPACT_OTHER_REPORT_ITEM_DELIM, -1);
        if(transValues.length != VAR_IMPACT_OTHER_REPORT_ITEM_COUNT)
            return null;

        return new AltTranscriptReportableInfo(
                transValues[0], transValues[1], transValues[2], transValues[3], CodingEffect.valueOf(transValues[4]));
    }
}
