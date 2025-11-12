package com.hartwig.hmftools.orange.report.finding;

import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.finding.Fusion;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FindingKeys {

    public static String findingKey(@NotNull PurpleVariant variant, @NotNull PurpleTranscriptImpact transcriptImpact) {
        return String.format("smallVariant[%s %s]", variant.gene(), impact(transcriptImpact));
    }

    public static String findingKey(PurpleGainDeletion gainDeletion) {
        return String.format("copyNumber[%s %s]", gainDeletion.gene(), gainDeletion.interpretation().name());
    }

    public static String findingKey(LinxBreakend breakend) {
        return String.format("disruption[%s]", breakend.gene());
    }

    public static String findingKey(LinxFusion fusion) {
        return String.format("fusion[%s %s]", fusion.geneStart(), fusion.geneEnd());
    }

    public static String findingKey(VirusInterpreterEntry virus) {
        String label = virus.interpretation() + (virus.interpretation() == VirusInterpretation.HPV ? " (" + virus.name() + ")" : "");
        return String.format("virus[%s]", label);
    }

    /*
    public static String findingKey(MicrosatelliteStability microsatelliteStability) {
        return MolecularCharacteristicEvents.MICROSATELLITE_UNSTABLE;
    }

    public static String findingKey(HomologousRecombination homologousRecombination) {
        return MolecularCharacteristicEvents.HOMOLOGOUS_RECOMBINATION_DEFICIENT;
    }

    public static String findingKey(TumorMutationalBurden tumorMutationalBurden) {
        return MolecularCharacteristicEvents.HIGH_TUMOR_MUTATIONAL_BURDEN;
    }

    public static String findingKey(TumorMutationalLoad tumorMutationalLoad) {
        return MolecularCharacteristicEvents.HIGH_TUMOR_MUTATIONAL_LOAD;
    }
    */

    public static String findingKey(@NotNull CuppaPrediction cuppaPrediction) {
        return String.format("predictedTumorOrigin[%s]", cuppaPrediction.cancerType());
    }

    private static String impact(@NotNull PurpleTranscriptImpact transcriptImpact) {
        return determineVariantAnnotation(transcriptImpact.hgvsCodingImpact(),
                transcriptImpact.hgvsProteinImpact(),
                transcriptImpact.effects().stream().map(Enum::toString).collect(Collectors.joining("&")),
                transcriptImpact.codingEffect() == PurpleCodingEffect.SPLICE,
                transcriptImpact.effects().contains(PurpleVariantEffect.UPSTREAM_GENE));
    }

    @NotNull
    public static String determineVariantAnnotation(@Nullable String hgvsCoding, @Nullable String hgvsProtein, @NotNull String effects,
            boolean isSplice, boolean isUpstream) {
        if (!Strings.isNullOrEmpty(hgvsProtein) && !hgvsProtein.equals("p.?")) {
            return hgvsProtein;
        }

        if (!Strings.isNullOrEmpty(hgvsCoding)) {
            return isSplice ? hgvsCoding + " splice" : hgvsCoding;
        }

        return isUpstream ? "upstream" : effects;
    }
}
