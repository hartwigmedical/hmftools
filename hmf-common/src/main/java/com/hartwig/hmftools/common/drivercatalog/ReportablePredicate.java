package com.hartwig.hmftools.common.drivercatalog;

import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class ReportablePredicate implements Predicate<SomaticVariant> {

    private final Map<String, DriverGene> driverGeneMap;

    public ReportablePredicate(@NotNull final DriverCategory type, @NotNull final DriverGenePanel driverGenePanel) {
        this.driverGeneMap = driverGenePanel.driverGenes()
                .stream()
                .filter(x -> x.likelihoodType().equals(type) && x.reportVariant())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    @Override
    public boolean test(final SomaticVariant variant) {
        final DriverGene driverGene = driverGeneMap.get(variant.gene());
        if (driverGene == null) {
            return false;
        }

        if (variant.isHotspot() && driverGene.reportHotspot()) {
            return true;
        }

        DriverImpact impact = DriverImpact.select(variant);
        switch (impact) {
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
