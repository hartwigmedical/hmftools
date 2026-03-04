package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.finding.datamodel.DriverSource;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

final class FindingKeys
{
    public static String smallVariant(DriverSource sampleType, PurpleVariant variant,
            PurpleTranscriptImpact transcriptImpact, boolean isCanonical)
    {
        return String.format("smallVariant[%s %s %s:%d %s %s]", sampleType, geneTranscriptLabel(variant.gene(), isCanonical,
                transcriptImpact.transcript()), variant.chromosome(), variant.position(), variant.ref(), variant.alt());
    }

    public static String gainDeletion(DriverSource sampleType, String gene, PurpleDriverType purpleDriverType,
            boolean isCanonical, String transcriptId)
    {
        return String.format("gainDeletion[%s %s %s]", sampleType, geneTranscriptLabel(gene, isCanonical, transcriptId),
                purpleDriverType.name());
    }

    public static String disruption(DriverSource sampleType, LinxBreakend breakend)
    {
        return String.format("disruption[%s %s %d]",
                sampleType,
                geneTranscriptLabel(breakend.gene(), breakend.isCanonical(), breakend.transcript()),
                breakend.svId());
    }

    public static String fusion(DriverSource sampleType, LinxFusion fusion)
    {
        return String.format("fusion[%s %s %s]", sampleType, fusion.geneStart(), fusion.geneEnd());
    }

    public static String virus(VirusInterpreterEntry virus)
    {
        String label = virus.interpretation() + (virus.interpretation() == VirusInterpretation.HPV ? " (" + virus.name() + ")" : "");
        return String.format("virus[%s]", label);
    }

    public static String hlaAllele(LilacAllele allele)
    {
        return String.format("hlaAllele[%s]", allele.allele());
    }

    public static String chromosomeArmCopyNumber(String chromosome, String arm)
    {
        return String.format("chrArmCopyNumber[%s%s]", chromosome, arm);
    }

    public static String microsatelliteStability(PurpleMicrosatelliteStatus status)
    {
        return String.format("microsatelliteStability[%s]", status.name());
    }

    public static String homologousRecombination(ChordStatus status)
    {
        return String.format("homologousRecombination[%s]", status.name());
    }

    public static String tumorMutationStatus(PurpleTumorMutationalStatus tmbStatus, PurpleTumorMutationalStatus tmlStatus)
    {
        return String.format("tumorMutationStatus[TMB_%s TML_%s]", tmbStatus, tmlStatus);
    }

    public static String predictedTumorOrigin(String cancerType)
    {
        return String.format("predictedTumorOrigin[%s]", cancerType);
    }

    public static String pharmacoGenotype(String gene, String allele)
    {
        return String.format("pharmacoGenotype[%s:%s]", gene, allele);
    }

    // only show transcript ID for non canonical transcripts
    private static String geneTranscriptLabel(String gene, boolean isCanonical, String transcriptId)
    {
        return isCanonical ? gene : String.format("%s(%s)", gene, transcriptId);
    }
}
