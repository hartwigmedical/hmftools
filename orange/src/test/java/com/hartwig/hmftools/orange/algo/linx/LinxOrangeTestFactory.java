package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.datamodel.driver.DriverInterpretation.HIGH;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LinxOrangeTestFactory
{
    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder()
    {
        return ImmutableLinxSvAnnotation.builder()
                .from(LinxConversion.convert(LinxTestFactory.svAnnotationBuilder().build()));
    }

    public static ImmutableLinxFusion.Builder fusionBuilder()
    {
        return ImmutableLinxFusion.builder()
                .geneUp("GENE_UP")
                .contextUp("")
                .transcriptUp("")
                .geneDown("GENE_DOWN")
                .contextDown("")
                .transcriptDown("")
                .reportedType(LinxFusionType.EXON_DEL_DUP)
                .phased(FusionPhasedType.INFRAME)
                .driverInterpretation(HIGH)
                .fusedExonUp(1)
                .fusedExonDown(2)
                .domainsKept("")
                .domainsLost("")
                .junctionCopyNumber(1)
                .chainLinks(0)
                .chainTerminated(false);
    }

    public static ImmutableLinxBreakend.Builder breakendBuilder()
    {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .gene("GENE")
                .chromosome("1")
                .chromosomeBand("2")
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .geneOrientation(LinxGeneOrientation.UPSTREAM)
                .disruptive(false)
                .reportedStatus(ReportedStatus.NONE)
                .undisruptedCopyNumber(0D)
                .type(LinxBreakendType.BND)
                .regionType(TranscriptRegionType.UNKNOWN)
                .codingType(TranscriptCodingType.UNKNOWN)
                .nextSpliceExonRank(0)
                .orientation(0)
                .exonUp(0)
                .exonDown(0)
                .junctionCopyNumber(0D)
                .driverType(LinxDriverType.DISRUPTION)
                .driverLikelihood(0);
    }

    public static ImmutableHomozygousDisruption.Builder homozygousDisruptionBuilder()
    {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false);
    }

    public static ImmutableLinxHomozygousDisruption.Builder linxHomozygousDisruptionBuilder()
    {
        return ImmutableLinxHomozygousDisruption.builder()
                .from(LinxConversion.convert(homozygousDisruptionBuilder().build()));
    }

}
