package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.datamodel.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverEventType;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

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
}
