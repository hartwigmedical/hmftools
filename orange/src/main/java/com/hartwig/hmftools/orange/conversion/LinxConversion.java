package com.hartwig.hmftools.orange.conversion;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.FusionReportableReason;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverEventType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
import com.hartwig.hmftools.orange.algo.linx.HomozygousDisruption;

public final class LinxConversion
{
    public static LinxDriver convert(final com.hartwig.hmftools.common.linx.LinxDriver driver)
    {
        return ImmutableLinxDriver.builder()
                .gene(driver.gene())
                .type(LinxDriverEventType.valueOf(driver.eventType().name()))
                .build();
    }

    public static LinxSvAnnotation convert(final com.hartwig.hmftools.common.linx.LinxSvAnnotation linxSvAnnotation)
    {
        return ImmutableLinxSvAnnotation.builder()
                .vcfId(linxSvAnnotation.vcfIdStart())
                .svId(linxSvAnnotation.svId())
                .clusterId(linxSvAnnotation.clusterId())
                .clusterReason(linxSvAnnotation.clusterReason())
                .fragileSiteStart(linxSvAnnotation.fragileSiteStart())
                .fragileSiteEnd(linxSvAnnotation.fragileSiteEnd())
                .isFoldback(linxSvAnnotation.isFoldback())
                .lineTypeStart(linxSvAnnotation.lineTypeStart())
                .lineTypeEnd(linxSvAnnotation.lineTypeEnd())
                .junctionCopyNumberMin(linxSvAnnotation.junctionCopyNumberMin())
                .junctionCopyNumberMax(linxSvAnnotation.junctionCopyNumberMax())
                .geneStart(linxSvAnnotation.geneStart())
                .geneEnd(linxSvAnnotation.geneEnd())
                .localTopologyIdStart(linxSvAnnotation.localTopologyIdStart())
                .localTopologyIdEnd(linxSvAnnotation.localTopologyIdEnd())
                .localTopologyStart(linxSvAnnotation.localTopologyStart())
                .localTopologyEnd(linxSvAnnotation.localTopologyEnd())
                .localTICountStart(linxSvAnnotation.localTICountStart())
                .localTICountEnd(linxSvAnnotation.localTICountEnd())
                .build();
    }

    public static LinxHomozygousDisruption convert(final DriverCatalog homozygousDisruption)
    {
        String typeStr = homozygousDisruption.driver() == DriverType.HOM_DEL_DISRUPTION ? "DEL" : "DUP";

        return ImmutableLinxHomozygousDisruption.builder()
                .gene(homozygousDisruption.gene())
                .chromosome(homozygousDisruption.chromosome())
                .chromosomeBand(homozygousDisruption.chromosomeBand())
                .transcript(homozygousDisruption.transcript())
                .isCanonical(homozygousDisruption.isCanonical())
                .driverInterpretation(DriverInterpretation.interpret(homozygousDisruption.driverLikelihood()))
                .type(typeStr)
                .build();
    }
}
