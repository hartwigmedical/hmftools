package com.hartwig.hmftools.orange.conversion;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.virus.VirusConstants;
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
import com.hartwig.hmftools.datamodel.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;

public final class OrangeConversion
{
    @NotNull
    public static com.hartwig.hmftools.datamodel.flagstat.Flagstat convert(Flagstat flagstat)
    {
        return ImmutableFlagstat.builder()
                .uniqueReadCount(flagstat.uniqueReadCount())
                .secondaryCount(flagstat.secondaryCount())
                .supplementaryCount(flagstat.supplementaryCount())
                .mappedProportion(flagstat.mappedProportion())
                .build();
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.metrics.WGSMetrics convert(WGSMetrics wgsMetrics)
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
    public static LilacRecord convert(LilacSummaryData lilacSummaryData)
    {
        return ImmutableLilacRecord.builder()
                .qc(lilacSummaryData.qc())
                .alleles(() -> lilacSummaryData.alleles().stream().map(OrangeConversion::convert).iterator())
                .build();
    }

    @NotNull
    public static LilacAllele convert(com.hartwig.hmftools.common.hla.LilacAllele allele)
    {
        return ImmutableLilacAllele.builder()
                .allele(allele.allele())
                .tumorCopyNumber(allele.tumorCopyNumber())
                .somaticMissense(allele.somaticMissense())
                .somaticNonsenseOrFrameshift(allele.somaticNonsenseOrFrameshift())
                .somaticSplice(allele.somaticSplice())
                .somaticSynonymous(allele.somaticSynonymous())
                .somaticInframeIndel(allele.somaticInframeIndel())
                .refFragments(allele.refFragments())
                .tumorFragments(allele.tumorFragments())
                .rnaFragments(allele.rnaFragments())
                .build();
    }

    @NotNull
    public static VirusInterpreterData convert(com.hartwig.hmftools.common.virus.VirusInterpreterData interpreterData)
    {
        return ImmutableVirusInterpreterData.builder()
                .allViruses(ConversionUtil.mapToIterable(interpreterData.allViruses(), OrangeConversion::convert))
                .reportableViruses(ConversionUtil.mapToIterable(interpreterData.reportableViruses(), OrangeConversion::convert))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static AnnotatedVirus convert(com.hartwig.hmftools.common.virus.AnnotatedVirus annotatedVirus)
    {
        VirusConstants interpretation = annotatedVirus.interpretation();
        return ImmutableAnnotatedVirus.builder()
                .name(annotatedVirus.name())
                .qcStatus(VirusBreakendQCStatus.valueOf(annotatedVirus.qcStatus().name()))
                .integrations(annotatedVirus.integrations())
                .interpretation(interpretation != null ? VirusInterpretation.valueOf(interpretation.toString()) : null)
                .percentageCovered(annotatedVirus.percentageCovered())
                .reported(annotatedVirus.reported())
                .meanCoverage(annotatedVirus.meanCoverage())
                .virusDriverLikelihoodType(VirusLikelihoodType.valueOf(annotatedVirus.virusDriverLikelihoodType().name()))
                .expectedClonalCoverage(annotatedVirus.expectedClonalCoverage())
                .build();
    }

    @NotNull
    public static ChordRecord convert(ChordData chordData)
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
    public static PeachGenotype convert(com.hartwig.hmftools.common.peach.PeachGenotype peachGenotype)
    {
        return ImmutablePeachGenotype.builder()
                .gene(peachGenotype.gene())
                .haplotype(peachGenotype.haplotype())
                .function(peachGenotype.function())
                .linkedDrugs(peachGenotype.linkedDrugs())
                .urlPrescriptionInfo(peachGenotype.urlPrescriptionInfo())
                .panelVersion(peachGenotype.panelVersion())
                .repoVersion(peachGenotype.repoVersion())
                .build();
    }

    @NotNull
    public static SignatureAllocation convert(com.hartwig.hmftools.common.sigs.SignatureAllocation signatureAllocation)
    {
        return ImmutableSignatureAllocation.builder()
                .signature(signatureAllocation.signature())
                .allocation(signatureAllocation.allocation())
                .percent(signatureAllocation.percent())
                .build();
    }
}
