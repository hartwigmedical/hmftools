package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.field.FieldInfo.findField;

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
import com.hartwig.hmftools.compar.common.field.FieldInfo;
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

    public PurityData(final PurityContext purityContext, final List<FieldInfo> fields)
    {
        mValues = Maps.newHashMap();

        addDoubleValue(PurityComparer.PurityFields.Purity, purityContext.bestFit().purity(), fields);
        addDoubleValue(PurityComparer.PurityFields.Ploidy, purityContext.bestFit().ploidy(), fields);
        addDoubleValue(PurityComparer.PurityFields.Contamination, purityContext.qc().contamination(), fields);
        addDoubleValue(PurityComparer.PurityFields.TmbPerMb, purityContext.tumorMutationalBurdenPerMb(), fields);
        addDoubleValue(PurityComparer.PurityFields.MsIndelsPerMb, purityContext.microsatelliteIndelsPerMb(), fields);

        addIntValue(PurityComparer.PurityFields.Tml, purityContext.tumorMutationalLoad(), fields);
        addIntValue(PurityComparer.PurityFields.CopyNumberSegments, purityContext.qc().copyNumberSegments(), fields);
        addIntValue(
                PurityComparer.PurityFields.UnsupportedCopyNumberSegments, purityContext.qc().unsupportedCopyNumberSegments(),
                fields);
        addIntValue(PurityComparer.PurityFields.SvTmb, purityContext.svTumorMutationalBurden(), fields);

        addStringValue(
                PurityComparer.PurityFields.QcStatus,
                purityContext.qc().status().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);
        addStringValue(PurityComparer.PurityFields.Gender, purityContext.gender().toString(), fields);
        addStringValue(
                PurityComparer.PurityFields.GermlineAberrations,
                purityContext.qc().germlineAberrations().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);
        addStringValue(PurityComparer.PurityFields.FitMethod, purityContext.method().toString(), fields);
        addStringValue(PurityComparer.PurityFields.MsStatus, purityContext.microsatelliteStatus().toString(), fields);
        addStringValue(PurityComparer.PurityFields.TmbStatus, purityContext.tumorMutationalBurdenStatus().toString(), fields);
        addStringValue(PurityComparer.PurityFields.TmlStatus, purityContext.tumorMutationalLoadStatus().toString(), fields);

        addDoubleValue(PurityComparer.PurityFields.TincLevel, purityContext.qc().tincLevel(), fields);

        Purity = purityContext;
    }

    public PurityData(final List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
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
                    addDoubleValue(field, Double.valueOf(truthsetValue.Value), fields);
                    break;

                case Tml:
                case CopyNumberSegments:
                case UnsupportedCopyNumberSegments:
                case SvTmb:
                    addIntValue(field, Integer.valueOf(truthsetValue.Value), fields);
                    break;

                case QcStatus:
                case Gender:
                case GermlineAberrations:
                case FitMethod:
                case MsStatus:
                case TmbStatus:
                case TmlStatus:
                    addStringValue(field, truthsetValue.Value, fields);
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

    // TODO: make generic, taking a string for field name
    private void addDoubleValue(final PurityComparer.PurityFields field, final Double value, final List<FieldInfo> fields)
    {
        mValues.put(field.name(), new DoubleFieldValue(findField(field.name(), fields), value));
    }

    private void addIntValue(final PurityComparer.PurityFields field, final Integer value, final List<FieldInfo> fields)
    {
        mValues.put(field.name(), new IntFieldValue(findField(field.name(), fields), value));
    }

    private void addStringValue(final PurityComparer.PurityFields field, final String value, final List<FieldInfo> fields)
    {
        mValues.put(field.name(), new StringFieldValue(findField(field.name(), fields), value));
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
