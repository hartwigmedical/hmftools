package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.chord.ImmutableChordRecord;
import com.hartwig.hmftools.datamodel.flagstat.ImmutableFlagstat;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.metrics.ImmutableWGSMetrics;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.datamodel.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;

public final class OrangeConversion
{
    public static com.hartwig.hmftools.datamodel.flagstat.Flagstat convert(final BamFlagStats flagstat)
    {
        return ImmutableFlagstat.builder()
                .uniqueReadCount(flagstat.uniqueReadCount())
                .secondaryCount(flagstat.secondaryCount())
                .supplementaryCount(flagstat.supplementaryCount())
                .mappedProportion(flagstat.mappedProportion())
                .build();
    }

    public static com.hartwig.hmftools.datamodel.metrics.WGSMetrics convert(final BamMetricSummary metricsSummary)
    {
        return ImmutableWGSMetrics.builder()
                .meanCoverage(metricsSummary.meanCoverage())
                .sdCoverage(metricsSummary.sdCoverage())
                .medianCoverage(metricsSummary.medianCoverage())
                .madCoverage(metricsSummary.madCoverage())
                .pctExcAdapter(0.0)
                .pctExcMapQ(metricsSummary.lowMapQualPercent())
                .pctExcDupe(metricsSummary.duplicatePercent())
                .pctExcUnpaired(metricsSummary.unmappedPercent())
                .pctExcBaseQ(metricsSummary.lowBaseQualPercent())
                .pctExcOverlap(metricsSummary.overlappingReadPercent())
                .pctExcCapped(metricsSummary.cappedCoveragePercent())
                .pctExcTotal(metricsSummary.totalFilteredPercent())
                .build();
    }

    public static OrangeDoidNode convert(final DoidNode doidNode)
    {
        return ImmutableOrangeDoidNode.builder().doid(doidNode.doid()).doidTerm(doidNode.doidTerm()).build();
    }

    public static LilacRecord convert(final LilacSummaryData lilacSummaryData, boolean hasRef, boolean hasRna)
    {
        return ImmutableLilacRecord.builder()
                .qc(lilacSummaryData.qc())
                .alleles(() -> lilacSummaryData.alleles()
                        .stream()
                        .map(allele -> OrangeConversion.convert(allele, hasRef, hasRna))
                        .iterator())
                .build();
    }

    public static LilacAllele convert(final com.hartwig.hmftools.common.hla.LilacAllele allele, boolean hasRef, boolean hasRna)
    {
        return ImmutableLilacAllele.builder()
                .allele(allele.allele())
                .tumorCopyNumber(allele.tumorCopyNumber())
                .refFragments(hasRef ? allele.refFragments() : null)
                .tumorFragments(allele.tumorFragments())
                .rnaFragments(hasRna ? allele.rnaFragments() : null)
                .somaticMissense(allele.somaticMissense())
                .somaticNonsenseOrFrameshift(allele.somaticNonsenseOrFrameshift())
                .somaticSplice(allele.somaticSplice())
                .somaticSynonymous(allele.somaticSynonymous())
                .somaticInframeIndel(allele.somaticInframeIndel())
                .build();
    }

    public static VirusInterpreterEntry convert(final com.hartwig.hmftools.common.virus.AnnotatedVirus annotatedVirus)
    {
        VirusType interpretation = annotatedVirus.interpretation();
        return ImmutableVirusInterpreterEntry.builder()
                .name(annotatedVirus.name())
                .qcStatus(VirusBreakendQCStatus.valueOf(annotatedVirus.qcStatus().name()))
                .integrations(annotatedVirus.integrations())
                .interpretation(interpretation != null ? VirusInterpretation.valueOf(interpretation.name()) : null)
                .percentageCovered(annotatedVirus.percentageCovered())
                .meanCoverage(annotatedVirus.meanCoverage())
                .expectedClonalCoverage(annotatedVirus.expectedClonalCoverage())
                .reported(annotatedVirus.reported())
                .driverLikelihood(VirusLikelihoodType.valueOf(annotatedVirus.virusDriverLikelihoodType().name()))
                .build();
    }

    public static ChordRecord convert(final ChordData chordData)
    {
        return ImmutableChordRecord.builder()
                .brca1Value(chordData.BRCA1Value())
                .brca2Value(chordData.BRCA2Value())
                .hrdValue(chordData.hrdValue())
                .hrStatus(ChordStatus.valueOf(chordData.hrStatus().name()))
                .hrdType(chordData.hrdType())
                .build();
    }

    public static PeachGenotype convert(final com.hartwig.hmftools.common.peach.PeachGenotype peachGenotype)
    {
        return ImmutablePeachGenotype.builder()
                .gene(peachGenotype.gene())
                .allele(peachGenotype.allele())
                .alleleCount(peachGenotype.alleleCount())
                .haplotype(peachGenotype.haplotype())
                .function(peachGenotype.function())
                .linkedDrugs(peachGenotype.linkedDrugs())
                .urlPrescriptionInfo(peachGenotype.urlPrescriptionInfo())
                .panelVersion(peachGenotype.panelVersion())
                .repoVersion(peachGenotype.repoVersion())
                .build();
    }
}
