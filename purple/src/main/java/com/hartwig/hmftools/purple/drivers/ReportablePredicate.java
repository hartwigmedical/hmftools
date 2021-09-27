package com.hartwig.hmftools.purple.drivers;

import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class ReportablePredicate implements Predicate<SomaticVariant>
{
    public static final int MAX_ONCO_REPEAT_COUNT = 7;

    private final int maxRepeatCount;
    private final Map<String, DriverGene> driverGeneMap;

    public ReportablePredicate(@NotNull final DriverCategory type, @NotNull final DriverGenePanel driverGenePanel)
    {
        this.driverGeneMap = driverGenePanel.driverGenes()
                .stream()
                .filter(x -> x.likelihoodType().equals(type) && x.reportSomatic())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));
        this.maxRepeatCount = type == DriverCategory.ONCO ? MAX_ONCO_REPEAT_COUNT : -1;
    }

    @Override
    public boolean test(final SomaticVariant variant)
    {
        final DriverGene driverGene = driverGeneMap.get(variant.gene());
        if(driverGene == null)
        {
            return false;
        }

        if(variant.type().equals(VariantType.INDEL) && maxRepeatCount > 0 && variant.repeatCount() > maxRepeatCount)
        {
            return false;
        }

        if(variant.isHotspot() && driverGene.reportSomaticHotspot())
        {
            return true;
        }

        DriverImpact impact = DriverImpact.select(variant);
        switch(impact)
        {
            case NONSENSE:
            case FRAMESHIFT:
                return driverGene.reportNonsenseAndFrameshift();
            case SPLICE:
                return driverGene.reportSplice();
            case MISSENSE:
            case INFRAME:
                return driverGene.reportMissenseAndInframe();
            default:
                return false;
        }
    }
}
