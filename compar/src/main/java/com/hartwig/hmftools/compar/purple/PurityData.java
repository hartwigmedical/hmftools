package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.DoubleFieldValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldValue;
import com.hartwig.hmftools.compar.common.field.IntFieldValue;
import com.hartwig.hmftools.compar.common.field.StringFieldValue;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    private final Map<String,FieldValue> mValues;

    @Deprecated
    public PurityData(final PurityContext purityContext)
    {
        mValues = Maps.newHashMap();
        Purity = purityContext;
    }

    public PurityData(final PurityContext purityContext, final Map<String,FieldCheck> fieldCheckMap)
    {
        mValues = Maps.newHashMap();

        mValues.put(PurityComparer.PurityFields.Purity.toString(),
                purity(PurityComparer.PurityFields.Purity, purityContext.bestFit().purity(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.Ploidy.toString(),
                ploidy(PurityComparer.PurityFields.Ploidy, purityContext.bestFit().ploidy(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.Contamination.toString(),
                contamination(PurityComparer.PurityFields.Contamination, purityContext.qc().contamination(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.TmbPerMb.toString(),
                tmbPerMb(PurityComparer.PurityFields.TmbPerMb, purityContext.tumorMutationalBurdenPerMb(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.MsIndelsPerMb.toString(),
                msiIndels(PurityComparer.PurityFields.MsIndelsPerMb, purityContext.microsatelliteIndelsPerMb(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.Tml.toString(),
                tml(PurityComparer.PurityFields.Tml, purityContext.tumorMutationalLoad(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.CopyNumberSegments.toString(),
                cnSegments(PurityComparer.PurityFields.CopyNumberSegments, purityContext.qc().copyNumberSegments(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.UnsupportedCopyNumberSegments.toString(),
                unsupportedCnSegs(PurityComparer.PurityFields.UnsupportedCopyNumberSegments,
                        purityContext.qc().unsupportedCopyNumberSegments(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.SvTmb.toString(),
                svTmb(PurityComparer.PurityFields.SvTmb, purityContext.svTumorMutationalBurden(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.QcStatus.toString(),
                stringField(
                        PurityComparer.PurityFields.QcStatus,
                        purityContext.qc().status().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                        fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.Gender.toString(),
                stringField(PurityComparer.PurityFields.Gender, purityContext.gender().toString(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.GermlineAberrations.toString(),
                stringField(
                        PurityComparer.PurityFields.GermlineAberrations,
                        purityContext.qc().germlineAberrations().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                        fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.FitMethod.toString(),
                stringField(PurityComparer.PurityFields.FitMethod, purityContext.method().toString(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.MsStatus.toString(),
                stringField(PurityComparer.PurityFields.MsStatus, purityContext.microsatelliteStatus().toString(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.TmbStatus.toString(),
                stringField(PurityComparer.PurityFields.TmbStatus, purityContext.tumorMutationalBurdenStatus().toString(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.TmlStatus.toString(),
                stringField(PurityComparer.PurityFields.TmlStatus, purityContext.tumorMutationalLoadStatus().toString(), fieldCheckMap));

        mValues.put(PurityComparer.PurityFields.TincLevel.toString(),
                tincLevel(PurityComparer.PurityFields.TincLevel, purityContext.qc().tincLevel(), fieldCheckMap));

        Purity = purityContext;
    }

    public PurityData(final List<TruthsetValue> truthsetValues, final Map<String,FieldCheck> fieldCheckMap)
    {
        mValues = Maps.newHashMap();

        for(TruthsetValue truthsetValue : truthsetValues)
        {
            // TODO: validate prior to creating fields
            String fieldName = truthsetValue.FieldName;
            PurityComparer.PurityFields field = PurityComparer.PurityFields.valueOf(fieldName);

            switch(field)
            {
                case Purity:
                    mValues.put(fieldName, purity(field, Double.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case Ploidy:
                    mValues.put(fieldName, ploidy(field, Double.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case Contamination:
                    mValues.put(fieldName, contamination(field, Double.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case TmbPerMb:
                    mValues.put(fieldName, tmbPerMb(field, Double.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case MsIndelsPerMb:
                    mValues.put(fieldName, msiIndels(field, Double.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case Tml:
                    mValues.put(fieldName, tml(field, Integer.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case CopyNumberSegments:
                    mValues.put(fieldName, cnSegments(field, Integer.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case UnsupportedCopyNumberSegments:
                    mValues.put(fieldName, unsupportedCnSegs(field, Integer.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case SvTmb:
                    mValues.put(fieldName, svTmb(field, Integer.valueOf(truthsetValue.Value), fieldCheckMap));
                    break;

                case QcStatus:
                case Gender:
                case GermlineAberrations:
                case FitMethod:
                case MsStatus:
                case TmbStatus:
                case TmlStatus:
                    mValues.put(fieldName, stringField(field, truthsetValue.Value, fieldCheckMap));
                    break;

                default:
                    CMP_LOGGER.error("unknown purity field: {}", truthsetValue);
            }
        }

        Purity = null;
    }

    @Override
    public Map<String,FieldValue> fieldValues() { return mValues; }

    @Override
    public List<String> fieldNames()
    {
        return Arrays.stream(PurityComparer.PurityFields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    private DoubleFieldValue purity(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.04, null),
                "%.2f");
    }

    private DoubleFieldValue ploidy(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.1, null),
                "%.2f");
    }

    private DoubleFieldValue contamination(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.005, null),
                "%.4f");
    }

    private DoubleFieldValue tmbPerMb(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.1, 0.05),
                "%.2f");
    }

    private DoubleFieldValue msiIndels(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.1, 0.05),
                "%.2f");
    }

    private IntFieldValue tml(final PurityComparer.PurityFields field, int value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new IntFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 1.0, 0.05));
    }

    private IntFieldValue cnSegments(final PurityComparer.PurityFields field, int value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new IntFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 5.0, 0.2));
    }

    private IntFieldValue unsupportedCnSegs(final PurityComparer.PurityFields field, int value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new IntFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 5.0, 0.2));
    }

    private IntFieldValue svTmb(final PurityComparer.PurityFields field, int value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new IntFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 5.0, 0.05));
    }

    private StringFieldValue stringField(final PurityComparer.PurityFields field, final String value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new StringFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), null, null));
    }

    private DoubleFieldValue tincLevel(final PurityComparer.PurityFields field, double value, final Map<String,FieldCheck> fieldCheckMap)
    {
        return new DoubleFieldValue(
                field.toString(), value, getOrMakeFieldCheck(fieldCheckMap, field.toString(), 0.1, null),
                "%.3f");
    }

    public CategoryType category() { return PURITY; }

    @Override
    public String key() { return ""; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public boolean supportTruthsetData() { return true; }
}
