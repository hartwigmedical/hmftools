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
import com.hartwig.hmftools.compar.common.field.ThresholdFieldCheck;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    private final Map<String, FieldValue> mValues;

    @Deprecated
    public PurityData(final PurityContext purityContext)
    {
        mValues = Maps.newHashMap();
        Purity = purityContext;
    }

    public PurityData(final PurityContext purityContext, final Map<String, FieldCheck> fieldCheckMap)
    {
        mValues = Maps.newHashMap();

        addDoubleValue(PurityComparer.PurityFields.Purity, purityContext.bestFit().purity(), fieldCheckMap);
        addDoubleValue(PurityComparer.PurityFields.Ploidy, purityContext.bestFit().ploidy(), fieldCheckMap);
        addDoubleValue(PurityComparer.PurityFields.Contamination, purityContext.qc().contamination(), fieldCheckMap);
        addDoubleValue(PurityComparer.PurityFields.TmbPerMb, purityContext.tumorMutationalBurdenPerMb(), fieldCheckMap);
        addDoubleValue(PurityComparer.PurityFields.MsIndelsPerMb, purityContext.microsatelliteIndelsPerMb(), fieldCheckMap);

        addIntValue(PurityComparer.PurityFields.Tml, purityContext.tumorMutationalLoad(), fieldCheckMap);
        addIntValue(PurityComparer.PurityFields.CopyNumberSegments, purityContext.qc().copyNumberSegments(), fieldCheckMap);
        addIntValue(
                PurityComparer.PurityFields.UnsupportedCopyNumberSegments, purityContext.qc().unsupportedCopyNumberSegments(),
                fieldCheckMap);
        addIntValue(PurityComparer.PurityFields.SvTmb, purityContext.svTumorMutationalBurden(), fieldCheckMap);

        addStringValue(
                PurityComparer.PurityFields.QcStatus,
                purityContext.qc().status().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fieldCheckMap);
        addStringValue(PurityComparer.PurityFields.Gender, purityContext.gender().toString(), fieldCheckMap);
        addStringValue(
                PurityComparer.PurityFields.GermlineAberrations,
                purityContext.qc().germlineAberrations().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fieldCheckMap);
        addStringValue(PurityComparer.PurityFields.FitMethod, purityContext.method().toString(), fieldCheckMap);
        addStringValue(PurityComparer.PurityFields.MsStatus, purityContext.microsatelliteStatus().toString(), fieldCheckMap);
        addStringValue(PurityComparer.PurityFields.TmbStatus, purityContext.tumorMutationalBurdenStatus().toString(), fieldCheckMap);
        addStringValue(PurityComparer.PurityFields.TmlStatus, purityContext.tumorMutationalLoadStatus().toString(), fieldCheckMap);

        addDoubleValue(PurityComparer.PurityFields.TincLevel, purityContext.qc().tincLevel(), fieldCheckMap);

        Purity = purityContext;
    }

    public PurityData(final List<TruthsetValue> truthsetValues, final Map<String, FieldCheck> fieldCheckMap)
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
                case Ploidy:
                case Contamination:
                case TmbPerMb:
                case MsIndelsPerMb:
                    addDoubleValue(field, Double.valueOf(truthsetValue.Value), fieldCheckMap);
                    break;

                case Tml:
                case CopyNumberSegments:
                case UnsupportedCopyNumberSegments:
                case SvTmb:
                    addIntValue(field, Integer.valueOf(truthsetValue.Value), fieldCheckMap);
                    break;

                case QcStatus:
                case Gender:
                case GermlineAberrations:
                case FitMethod:
                case MsStatus:
                case TmbStatus:
                case TmlStatus:
                    addStringValue(field, truthsetValue.Value, fieldCheckMap);
                    break;

                default:
                    CMP_LOGGER.error("unknown purity field: {}", truthsetValue);
            }
        }

        Purity = null;
    }

    @Override
    public Map<String, FieldValue> fieldValues() { return mValues; }

    @Override
    public List<String> fieldNames()
    {
        return Arrays.stream(PurityComparer.PurityFields.values()).map(x -> x.name()).collect(Collectors.toList());
    }

    private void addDoubleValue(final PurityComparer.PurityFields field, final Double value,
            final Map<String, FieldCheck> fieldCheckMap)
    {
        mValues.put(
                field.name(),
                new DoubleFieldValue(field.name(), value, determineThresholdFieldCheck(field, fieldCheckMap), field.FormatString));
    }

    private void addIntValue(final PurityComparer.PurityFields field, final Integer value,
            final Map<String, FieldCheck> fieldCheckMap)
    {
        mValues.put(field.name(), new IntFieldValue(field.name(), value, determineThresholdFieldCheck(field, fieldCheckMap)));
    }

    private void addStringValue(final PurityComparer.PurityFields field, final String value, final Map<String, FieldCheck> fieldCheckMap)
    {
        mValues.put(field.name(), new StringFieldValue(field.name(), value, determineFieldCheck(field, fieldCheckMap)));
    }

    private static ThresholdFieldCheck determineThresholdFieldCheck(final PurityComparer.PurityFields field,
            final Map<String, FieldCheck> fieldCheckMap)
    {
        return getOrMakeFieldCheck(fieldCheckMap, field.name(), field.DefaultAbsoluteThreshold, field.DefaultPercentageThreshold);
    }

    private static FieldCheck determineFieldCheck(final PurityComparer.PurityFields field,
            final Map<String, FieldCheck> fieldCheckMap)
    {
        FieldCheck override = fieldCheckMap.get(field.name());
        return override != null ? override : new FieldCheck(true);
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
