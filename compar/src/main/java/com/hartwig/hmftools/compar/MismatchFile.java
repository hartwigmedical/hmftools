package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_SAMPLE_ID;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.compar.common.CurationCache.Entry;
import com.hartwig.hmftools.compar.common.FieldDisplayInfo;
import com.hartwig.hmftools.compar.common.Mismatch;

public class MismatchFile
{
    public static final String FLD_MISMATCH_TYPE = "MismatchType";
    public static final String FLD_CURATION_STATUS = "CurationStatus";
    private static final String FLD_OLD_VALUES = "OldValues";
    private static final String FLD_NEW_VALUES = "NewValues";
    public static final String FLD_DIFFERENCES = "Differences";
    public static final String FLD_OLD_SOURCE = "OldSource";
    public static final String FLD_NEW_SOURCE = "NewSource";
    public static final String FLD_CURATION_TYPE = "CurationType";
    public static final String FLD_VERDICT = "Verdict";
    public static final String FLD_CURATION_COMMENT = "CurationComment";
    public static final String FLD_CURATED_BY = "CuratedBy";
    public static final String FLD_CURATED_DATE = "CuratedDate";

    private static final String FLD_PREV_COMPARISON = "PrevComparison";
    private static final String FLD_PREV_MISMATCH_TYPE = "PrevMismatchType";
    private static final String FLD_PREV_DIFFERENCES = "PrevDifferences";
    private static final String FLD_PREV_CURATION_TYPE = "PrevCurationType";
    private static final String FLD_PREV_VERDICT = "PrevVerdict";
    private static final String FLD_PREV_CURATION_COMMENT = "PrevCurationComment";
    private static final String FLD_PREV_CURATED_BY = "PrevCuratedBy";
    private static final String FLD_PREV_CURATED_DATE = "PrevCuratedDate";

    private static final String SOURCE_PAIR_DELIM = ">";

    // differences is list of the form: field(oldValue/newValue)
    public static String commonHeader(boolean includeSampleId, boolean includeCatagory, boolean includeCurationStatus)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(includeSampleId)
            sj.add(FLD_SAMPLE_ID.toString());

        if(includeCatagory)
            sj.add(FLD_CATEGORY);

        sj.add(FLD_MISMATCH_TYPE);

        if(includeCurationStatus)
            sj.add(FLD_CURATION_STATUS);

        sj.add(FLD_ITEM_KEY).add(FLD_DIFFERENCES);
        return sj.toString();
    }

    public static String header(boolean includeSampleId, boolean includeSourceLabels, boolean includeHistory)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(commonHeader(includeSampleId, true, includeSourceLabels || includeHistory));
        sj.add(FLD_OLD_VALUES);
        sj.add(FLD_NEW_VALUES);

        if(includeSourceLabels)
            sj.add(FLD_OLD_SOURCE).add(FLD_NEW_SOURCE);

        if(includeHistory)
        {
            sj.add(FLD_PREV_COMPARISON).add(FLD_PREV_MISMATCH_TYPE).add(FLD_PREV_DIFFERENCES).add(FLD_PREV_CURATION_TYPE);
            sj.add(FLD_PREV_VERDICT).add(FLD_PREV_CURATION_COMMENT).add(FLD_PREV_CURATED_BY).add(FLD_PREV_CURATED_DATE);
        }

        return sj.toString();
    }

    public static String commonTsv(boolean writeCategory, final Mismatch mismatch, final String curationStatus)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(writeCategory)
        {
            if(mismatch.OldItem != null)
                sj.add(mismatch.OldItem.category().toString());
            else
                sj.add(mismatch.NewItem.category().toString());
        }

        sj.add(mismatch.Type.toString());

        if(curationStatus != null)
            sj.add(curationStatus);

        if(mismatch.OldItem != null)
            sj.add(mismatch.OldItem.key());
        else
            sj.add(mismatch.NewItem.key());

        return sj.toString();
    }

    public static String toTsv(
            final Mismatch mismatch, boolean writeFieldValues, final List<FieldDisplayInfo> fieldDisplayValues, final String curationStatus)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(commonTsv(!writeFieldValues, mismatch, curationStatus));

        sj.add(String.join(ITEM_DELIM, mismatch.DiffValues));

        if(writeFieldValues)
        {
            for(FieldDisplayInfo fieldInfo : fieldDisplayValues)
            {
                sj.add(fieldInfo.oldValue() != null ? fieldInfo.oldValue() : "");
                sj.add(fieldInfo.newValue() != null ? fieldInfo.newValue() : "");
            }
        }
        else
        {
            sj.add(itemValues(mismatch.OldItem, fieldDisplayValues, true));
            sj.add(itemValues(mismatch.NewItem, fieldDisplayValues, false));
        }

        return sj.toString();
    }

    public static String historyTsv(final Entry previous)
    {
        if(previous == null)
            return TSV_DELIM.repeat(7);

        return String.join(TSV_DELIM, previous.oldSource() + SOURCE_PAIR_DELIM + previous.newSource(), previous.type().toString(),
                previous.differences(), previous.curation().toString(), previous.verdict().toString(), previous.comment(),
                previous.curatedBy(), previous.curatedDate());
    }

    private static String itemValues(final ComparableItem item, final List<FieldDisplayInfo> fieldDisplayValues, boolean useOld)
    {
        if(item == null)
            return "";

        StringJoiner displaySj = new StringJoiner(ITEM_DELIM);

        for(FieldDisplayInfo fieldInfo : fieldDisplayValues)
        {
            String value = useOld ? fieldInfo.oldValue() : fieldInfo.newValue();

            if(value != null) // supress values not present (ie in the truthset)
                displaySj.add(format("%s=%s", fieldInfo.name(), value));
        }

        return displaySj.toString();
    }
}
