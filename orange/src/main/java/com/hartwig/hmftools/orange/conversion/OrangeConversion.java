package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
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
    @NotNull
    public static com.hartwig.hmftools.datamodel.flagstat.Flagstat convert(@NotNull Flagstat flagstat)
    {
        return ImmutableFlagstat.builder()
                .uniqueReadCount(flagstat.uniqueReadCount())
                .secondaryCount(flagstat.secondaryCount())
                .supplementaryCount(flagstat.supplementaryCount())
                .mappedProportion(flagstat.mappedProportion())
                .build();
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.metrics.WGSMetrics convert(@NotNull WGSMetrics wgsMetrics)
    {
        return ImmutableWGSMetrics.builder()
                .meanCoverage(wgsMetrics.meanCoverage())
                .sdCoverage(wgsMetrics.sdCoverage())
                .medianCoverage(wgsMetrics.medianCoverage())
                .madCoverage(wgsMetrics.madCoverage())
                .pctExcAdapter(wgsMetrics.pctExcAdapter())
                .pctExcMapQ(wgsMetrics.pctExcMapQ())
                .pctExcDupe(wgsMetrics.pctExcDupe())
                .pctExcUnpaired(wgsMetrics.pctExcUnpaired())
                .pctExcBaseQ(wgsMetrics.pctExcBaseQ())
                .pctExcOverlap(wgsMetrics.pctExcOverlap())
                .pctExcCapped(wgsMetrics.pctExcCapped())
                .pctExcTotal(wgsMetrics.pctExcTotal())
                .build();
    }

    @NotNull
    public static OrangeDoidNode convert(@NotNull DoidNode doidNode)
    {
        return ImmutableOrangeDoidNode.builder().doid(doidNode.doid()).doidTerm(doidNode.doidTerm()).build();
    }

    @NotNull
    public static LilacRecord convert(@NotNull LilacSummaryData lilacSummaryData, boolean hasRef, boolean hasRna)
    {
        return ImmutableLilacRecord.builder()
                .qc(lilacSummaryData.qc())
                .alleles(() -> lilacSummaryData.alleles()
                        .stream()
                        .map(allele -> OrangeConversion.convert(allele, hasRef, hasRna))
                        .iterator())
                .build();
    }

    @NotNull
    public static LilacAllele convert(@NotNull com.hartwig.hmftools.common.hla.LilacAllele allele, boolean hasRef, boolean hasRna)
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

    @NotNull
    public static VirusInterpreterEntry convert(@NotNull com.hartwig.hmftools.common.virus.AnnotatedVirus annotatedVirus)
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

    @NotNull
    public static ChordRecord convert(@NotNull ChordData chordData)
    {
        return ImmutableChordRecord.builder()
                .brca1Value(chordData.BRCA1Value())
                .brca2Value(chordData.BRCA2Value())
                .hrdValue(chordData.hrdValue())
                .hrStatus(ChordStatus.valueOf(chordData.hrStatus().name()))
                .hrdType(chordData.hrdType())
                .build();
    }

    @NotNull
    public static PeachGenotype convert(@NotNull com.hartwig.hmftools.common.peach.PeachGenotype peachGenotype)
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
