package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class VirusData implements ComparableItem
{
    public final AnnotatedVirus Virus;

    protected static final String FLD_INTEGRATIONS = "Integrations";
    protected static final String FLD_MEAN_COVERAGE = "MeanCoverage";
    protected static final String FLD_DRIVER_LIKELIHOOD = "DriverLikelihood";

    VirusData(final AnnotatedVirus virus)
    {
        Virus = virus;
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
