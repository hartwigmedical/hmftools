package com.hartwig.hmftools.common.drivercatalog;

import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverLikelihoodType;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class OncoPredicate implements Predicate<SomaticVariant> {

    static final int MAX_REPEAT_COUNT = 7; //TODO: INVESTIGATE THIS
    private final Map<String, DriverGene> driverGeneMap;

    public OncoPredicate(@NotNull final DriverGenePanel driverGenePanel) {
        this.driverGeneMap = driverGenePanel.driverGenes()
                .stream()
                .filter(x -> x.likelihoodType().equals(DriverLikelihoodType.ONCO))
                .collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    @Override
    public boolean test(final SomaticVariant variant) {
        final DriverGene driverGene = driverGeneMap.get(variant.gene());
        if (driverGene == null) {
            return false;
        }

        //TODO: DOUBLE CHECK THIS
        if (variant.isHotspot()) { //) && driverGene.reportPromoterHotspots()) {
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
