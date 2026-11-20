package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.MUTATION_BURDEN;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class MutationBurdenData extends ComparableItem
{
    public final PurityContext Purity;

    public MutationBurdenData(final PurityContext purityContext, final List<FieldInfo> fields)
    {
        addDoubleValue(MutationBurdenComparer.Fields.TmbPerMb.toString(), purityContext.tumorMutationalBurdenPerMb(), fields);
        addDoubleValue(MutationBurdenComparer.Fields.MsIndelsPerMb.toString(), purityContext.microsatelliteIndelsPerMb(), fields);

        addIntValue(MutationBurdenComparer.Fields.Tml.toString(), purityContext.tumorMutationalLoad(), fields);

        addStringValue(MutationBurdenComparer.Fields.MsStatus.toString(), purityContext.microsatelliteStatus().toString(), fields);
        addStringValue(MutationBurdenComparer.Fields.TmbStatus.toString(), purityContext.tumorMutationalBurdenStatus().toString(), fields);
        addStringValue(MutationBurdenComparer.Fields.TmlStatus.toString(), purityContext.tumorMutationalLoadStatus().toString(), fields);

        Purity = purityContext;
    }

    public MutationBurdenData(final List<ExpectedValue> expectedValues, final List<FieldInfo> fields)
    {
        for(ExpectedValue expectedValue : expectedValues)
        {
            String fieldName = expectedValue.FieldName;
            MutationBurdenComparer.Fields field = MutationBurdenComparer.Fields.valueOf(fieldName);

            switch(field)
            {
                case TmbPerMb:
                case MsIndelsPerMb:
                    addDoubleValue(field.toString(), Double.valueOf(expectedValue.Value), fields);
                    break;

                case Tml:
                    addIntValue(field.toString(), Integer.valueOf(expectedValue.Value), fields);
                    break;

                case MsStatus:
                case TmbStatus:
                case TmlStatus:
                    addStringValue(field.toString(), expectedValue.Value, fields);
                    break;

                default:
                    CMP_LOGGER.error("unknown mutation burden field: {}", expectedValue);
            }
        }

        Purity = null;
    }

    public CategoryType category() { return MUTATION_BURDEN; }

    @Override
    public String key() { return ""; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        return true;
    }
}
