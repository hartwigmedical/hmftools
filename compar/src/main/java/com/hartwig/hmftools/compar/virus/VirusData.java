package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;

import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class VirusData extends ComparableItem
{
    public final AnnotatedVirus Virus;

    VirusData(final AnnotatedVirus virus, final List<FieldInfo> fields)
    {
        Virus = virus;

        addBoolValue(VirusComparer.Fields.Reported.toString(), virus.reported(), fields);
        addIntValue(VirusComparer.Fields.Integrations.toString(), virus.integrations(), fields);
        addDoubleValue(VirusComparer.Fields.MeanCoverage.toString(), virus.meanCoverage(), fields);
        addStringValue(VirusComparer.Fields.DriverLikelihood.toString(), virus.virusDriverLikelihoodType().toString(), fields);
    }

    @Override
    public CategoryType category() {
        return VIRUS;
    }

    @Override
    public String key()
    {
        return Virus.name();
    }

    @Override
    public boolean reportable()
    {
        return Virus.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final VirusData otherData = (VirusData) other;
        return Virus.name().equals(otherData.Virus.name());
    }
}
