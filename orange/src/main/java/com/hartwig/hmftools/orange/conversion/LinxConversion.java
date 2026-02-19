package com.hartwig.hmftools.orange.conversion;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.linx.FusionReportableReason;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
import com.hartwig.hmftools.orange.algo.linx.HomozygousDisruption;

import org.jetbrains.annotations.NotNull;

public final class LinxConversion
{
    public static LinxDriver convert(final com.hartwig.hmftools.common.linx.LinxDriver driver)
    {
        return ImmutableLinxDriver.builder()
                .gene(driver.gene())
                .type(LinxDriverType.valueOf(driver.eventType().name()))
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

    public static LinxFusion convert(final com.hartwig.hmftools.common.linx.LinxFusion linxFusion)
    {
        return ImmutableLinxFusion.builder()
                .geneStart(linxFusion.geneStart())
                .geneContextStart(linxFusion.geneContextStart())
                .geneTranscriptStart(linxFusion.geneTranscriptStart())
                .geneEnd(linxFusion.geneEnd())
                .geneContextEnd(linxFusion.geneContextEnd())
                .geneTranscriptEnd(linxFusion.geneTranscriptEnd())
                .reported(linxFusion.reported())
                .reportedType(LinxFusionType.valueOf(linxFusion.reportedType()))
                .unreportedReasons(convertUnreportableReasons(linxFusion.reportableReasons()))
                .phased(FusionPhasedType.valueOf(linxFusion.phased().name()))
                .driverLikelihood(FusionLikelihoodType.valueOf(linxFusion.likelihood().name()))
                .fusedExonUp(linxFusion.fusedExonUp())
                .fusedExonDown(linxFusion.fusedExonDown())
                .chainLinks(linxFusion.chainLinks())
                .chainTerminated(linxFusion.chainTerminated())
                .domainsKept(linxFusion.domainsKept())
                .domainsLost(linxFusion.domainsLost())
                .junctionCopyNumber(linxFusion.junctionCopyNumber())
                .build();
    }

    private static List<LinxUnreportableReason> convertUnreportableReasons(final List<FusionReportableReason> reasons)
    {
        return reasons.stream()
                .map(item -> switch(item)
                {
                    case OK -> LinxUnreportableReason.NONE;
                    case PROTEIN_DOMAINS -> LinxUnreportableReason.DISRUPTED_PROTEIN_DOMAINS;
                    case NMD -> LinxUnreportableReason.NONSENSE_MEDIATED_DECAY;
                    default -> LinxUnreportableReason.valueOf(item.name());
                }).collect(Collectors.toList());
    }

    public static LinxHomozygousDisruption convert(final HomozygousDisruption homozygousDisruption)
    {
        return ImmutableLinxHomozygousDisruption.builder()
                .gene(homozygousDisruption.gene())
                .chromosome(homozygousDisruption.chromosome())
                .chromosomeBand(homozygousDisruption.chromosomeBand())
                .transcript(homozygousDisruption.transcript())
                .isCanonical(homozygousDisruption.isCanonical())
                .build();
    }
}
