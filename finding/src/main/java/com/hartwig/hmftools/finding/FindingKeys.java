package com.hartwig.hmftools.finding;

import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FindingKeys
{
    public static String smallVariant(@NotNull DriverSource sampleType, @NotNull PurpleVariant variant,
            @NotNull PurpleTranscriptImpact transcriptImpact, boolean isCanonical)
    {
        return String.format("smallVariant[%s %s %s]", sampleType, geneTranscriptLabel(variant.gene(), isCanonical,
                transcriptImpact.transcript()), impact(transcriptImpact));
    }

    public static String gainDeletion(@NotNull DriverSource sampleType, String gene, CopyNumberInterpretation copyNumberInterpretation,
            boolean isCanonical, String transcriptId)
    {
        return String.format("gainDeletion[%s %s %s]", sampleType, geneTranscriptLabel(gene, isCanonical, transcriptId),
                copyNumberInterpretation.name());
    }

    public static String disruption(@NotNull DriverSource sampleType, LinxBreakend breakend)
    {
        return String.format("disruption[%s %s %d]",
                sampleType,
                geneTranscriptLabel(breakend.gene(), breakend.isCanonical(), breakend.transcript()),
                breakend.svId());
    }

    public static String fusion(@NotNull DriverSource sampleType, LinxFusion fusion)
    {
        return String.format("fusion[%s %s %s]", sampleType, fusion.geneStart(), fusion.geneEnd());
    }

    public static String virus(VirusInterpreterEntry virus)
    {
        String label = virus.interpretation() + (virus.interpretation() == VirusInterpretation.HPV ? " (" + virus.name() + ")" : "");
        return String.format("virus[%s]", label);
    }

    public static String hlaAllele(@NotNull LilacAllele allele)
    {
        return String.format("hlaAllele[%s]", allele.allele());
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

    public static String predictedTumorOrigin(@NotNull String cancerType)
    {
        return String.format("predictedTumorOrigin[%s]", cancerType);
    }

    public static String pharmacoGenotype(@NotNull String gene, @NotNull String allele)
    {
        return String.format("pharmacoGenotype[%s:%s]", gene, allele);
    }

    // only show transcript ID for non canonical transcripts
    private static String geneTranscriptLabel(String gene, boolean isCanonical, String transcriptId)
    {
        return isCanonical ? gene : String.format("%s(%s)", gene, transcriptId);
    }

    private static String impact(@NotNull PurpleTranscriptImpact transcriptImpact)
    {
        return determineVariantAnnotation(transcriptImpact.hgvsCodingImpact(),
                transcriptImpact.hgvsProteinImpact(),
                transcriptImpact.effects().stream().map(Enum::toString).collect(Collectors.joining("&")),
                transcriptImpact.codingEffect() == PurpleCodingEffect.SPLICE,
                transcriptImpact.effects().contains(PurpleVariantEffect.UPSTREAM_GENE));
    }

    @NotNull
    public static String determineVariantAnnotation(@Nullable String hgvsCoding, @Nullable String hgvsProtein, @NotNull String effects,
            boolean isSplice, boolean isUpstream)
    {
        if(hgvsProtein != null && !hgvsProtein.isEmpty() && !hgvsProtein.equals("p.?"))
        {
            return hgvsProtein;
        }

        if(hgvsCoding != null && !hgvsCoding.isEmpty())
        {
            return isSplice ? hgvsCoding + " splice" : hgvsCoding;
        }

        return isUpstream ? "upstream" : effects;
    }
}
