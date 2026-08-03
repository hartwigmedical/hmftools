package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class PurityData extends ComparableItem
{
    public final PurityContext Purity;

    public PurityData(final PurityContext purityContext, final List<FieldInfo> fields)
    {
        addDoubleValue(PurityComparer.Fields.Purity.toString(), purityContext.bestFit().purity(), fields);
        addDoubleValue(PurityComparer.Fields.Ploidy.toString(), purityContext.bestFit().ploidy(), fields);
        addDoubleValue(PurityComparer.Fields.Contamination.toString(), purityContext.qc().contamination(), fields);
        addDoubleValue(PurityComparer.Fields.TmbPerMb.toString(), purityContext.tumorMutationalBurdenPerMb(), fields);
        addDoubleValue(PurityComparer.Fields.MsIndelsPerMb.toString(), purityContext.microsatelliteIndelsPerMb(), fields);

        addIntValue(PurityComparer.Fields.Tml.toString(), purityContext.tumorMutationalLoad(), fields);
        addIntValue(PurityComparer.Fields.CopyNumberSegments.toString(), purityContext.qc().copyNumberSegments(), fields);

        addIntValue(
                PurityComparer.Fields.UnsupportedCopyNumberSegments.toString(), purityContext.qc().unsupportedCopyNumberSegments(),
                fields);

        addIntValue(PurityComparer.Fields.SvTmb.toString(), purityContext.svTumorMutationalBurden(), fields);

        addStringValue(
                PurityComparer.Fields.QcStatus.toString(),
                purityContext.qc().status().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);

        addStringValue(PurityComparer.Fields.Gender.toString(), purityContext.gender().toString(), fields);

        addStringValue(
                PurityComparer.Fields.GermlineAberrations.toString(),
                purityContext.qc().germlineAberrations().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);

        addStringValue(PurityComparer.Fields.FitMethod.toString(), purityContext.method().toString(), fields);
        addStringValue(PurityComparer.Fields.MsStatus.toString(), purityContext.microsatelliteStatus().toString(), fields);
        addStringValue(PurityComparer.Fields.TmbStatus.toString(), purityContext.tumorMutationalBurdenStatus().toString(), fields);
        addStringValue(PurityComparer.Fields.TmlStatus.toString(), purityContext.tumorMutationalLoadStatus().toString(), fields);

        addDoubleValue(PurityComparer.Fields.TincLevel.toString(), purityContext.qc().tincLevel(), fields);

        Purity = purityContext;
    }

    public PurityData(final List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        for(TruthsetValue truthsetValue : truthsetValues)
        {
            // TODO: validate prior to creating fields
            String fieldName = truthsetValue.FieldName;
            PurityComparer.Fields field = PurityComparer.Fields.valueOf(fieldName);

            switch(field)
            {
                case Purity:
                case Ploidy:
                case Contamination:
                case TmbPerMb:
                case MsIndelsPerMb:
                    addDoubleValue(field.toString(), Double.valueOf(truthsetValue.Value), fields);
                    break;

                case Tml:
                case CopyNumberSegments:
                case UnsupportedCopyNumberSegments:
                case SvTmb:
                    addIntValue(field.toString(), Integer.valueOf(truthsetValue.Value), fields);
                    break;

                case QcStatus:
                case Gender:
                case GermlineAberrations:
                case FitMethod:
                case MsStatus:
                case TmbStatus:
                case TmlStatus:
                    addStringValue(field.toString(), truthsetValue.Value, fields);
                    break;

                default:
                    CMP_LOGGER.error("unknown purity field: {}", truthsetValue);
            }
        }

        Purity = null;
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
}
