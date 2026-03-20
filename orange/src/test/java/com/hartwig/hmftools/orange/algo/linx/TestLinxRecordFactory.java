package com.hartwig.hmftools.orange.algo.linx;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;

public final class TestLinxRecordFactory
{
    public static ImmutableLinxBreakend.Builder linxRecordBreakendBuilder()
    {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(1)
                .gene("GENE")
                .chromosome("1")
                .position(1000)
                .chromosomeBand("1P23")
                .transcript("TRANS")
                .isCanonical(true)
                .geneOrientation(LinxGeneOrientation.UPSTREAM)
                .disruptive(true)
                .reportedStatus(ReportedStatus.REPORTED)
                .undisruptedCopyNumber(1)
                .type(LinxBreakendType.BND)
                .regionType(TranscriptRegionType.INTRONIC)
                .codingType(TranscriptCodingType.CODING)
                .nextSpliceExonRank(2)
                .orientation(1)
                .exonUp(2)
                .exonDown(2)
                .junctionCopyNumber(1)
                .driverType(LinxDriverType.DISRUPTION)
                .driverLikelihood(0)
                .plotFilename(null);
    }

    public static ImmutableLinxFusion.Builder linxRecordFusionBuilder()
    {
        return ImmutableLinxFusion.builder()
                .geneUp("GENE01")
                .contextUp("")
                .transcriptUp("")
                .geneDown("GENE02")
                .contextDown("")
                .transcriptDown("")
                .reportedType(LinxFusionType.KNOWN_PAIR)
                .phased(FusionPhasedType.INFRAME)
                .driverInterpretation(DriverInterpretation.HIGH)
                .fusedExonUp(1)
                .fusedExonDown(2)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept("")
                .domainsLost("")
                .junctionCopyNumber(1)
                .rnaSupport(null)
                .plotFilename(null);
    }

    public static ImmutableLinxRecord.Builder linxRecordBuilder()
    {
        return ImmutableLinxRecord.builder()
                .somaticBreakends(Lists.newArrayList())
                .germlineBreakends(Lists.newArrayList())
                .fusions(Lists.newArrayList());
    }

    public static LinxRecord createMinimalTestLinxData()
    {
        return linxRecordBuilder().build();
    }
}
