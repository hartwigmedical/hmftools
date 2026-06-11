package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombination;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;

public final class FindingKeys
{
    public static String smallVariant(DriverSource driverSource, PurpleVariant variant,
            PurpleTranscriptImpact transcriptImpact, boolean isCanonical)
    {
        return String.format("smallVariant[%s %s %s:%d %s %s]", driverSource, geneTranscriptLabel(variant.gene(), isCanonical,
                transcriptImpact.transcript()), variant.chromosome(), variant.position(), variant.ref(), variant.alt());
    }

    public static String gainDeletion(DriverSource driverSource, String gene, PurpleDriverType driverType,
            boolean isCanonical, String transcriptId)
    {
        return String.format("gainDeletion[%s %s %s]", driverSource, geneTranscriptLabel(gene, isCanonical, transcriptId),
                driverType.name());
    }

    public static String disruption(DriverSource driverSource, LinxBreakend breakend)
    {
        return String.format("disruption[%s %s %d]",
                driverSource,
                geneTranscriptLabel(breakend.gene(), breakend.isCanonical(), breakend.transcript()),
                breakend.svId());
    }

    public static String fusion(DriverSource driverSource, LinxFusion fusion)
    {
        return String.format("fusion[%s %s %s]", driverSource, fusion.geneUp(), fusion.geneDown());
    }

    public static String virus(VirusInterpreterEntry virus)
    {
        String label = virus.interpretation() + (virus.interpretation() == VirusInterpretation.HPV ? " (" + virus.name() + ")" : "");
        return String.format("virus[%s]", label);
    }

    public static String chromosomeArmCopyNumber(String chromosome, String arm)
    {
        return String.format("chrArmCopyNumber[%s%s]", chromosome, arm);
    }

    public static String hlaAllele(LilacAllele allele, int suffix)
    {
        return String.format("hlaAllele[%s %d]", allele.allele(), suffix);
    }

    public static String microsatelliteStability(PurpleMicrosatelliteStatus status)
    {
        return String.format("microsatelliteStability[%s]", status.name());
    }

    public static String homologousRecombination(HomologousRecombination.Status status, HomologousRecombination.HrdCancerType cancerType)
    {
        return String.format("homologousRecombination[%s %s]", status.name(), cancerType.name());
    }

    public static String tumorMutationLoadStatus(PurpleTumorMutationalStatus status)
    {
        return String.format("tumorMutationLoadStatus[TML_%s]", status);
    }

    public static String tumorMutationBurdenStatus(PurpleTumorMutationalStatus status)
    {
        return String.format("tumorMutationBurdenStatus[TMB_%s]", status);
    }

    public static String predictedTumorOrigin(String cancerType)
    {
        return String.format("predictedTumorOrigin[%s]", cancerType);
    }

    public static String pharmacoGenotype(String gene, String allele)
    {
        return String.format("pharmacoGenotype[%s:%s]", gene, allele);
    }

    public static String rnaGeneExpression(String expressionType, GeneExpression geneExpression)
    {
        return String.format("rnaGeneExpression[%s %s]", expressionType, geneExpression.gene());
    }

    public static String rnaFusion(RnaFusion fusion)
    {
        return String.format("rnaFusion[%s %s:%d %s %s:%d]",
                fusion.geneStart(),
                fusion.chromosomeStart(),
                fusion.positionStart(),
                fusion.geneEnd(),
                fusion.chromosomeEnd(),
                fusion.positionEnd());
    }

    public static String novelSpliceJunction(NovelSpliceJunction spliceJunction)
    {
        return String.format("novelSpliceJunction[%s %s:%d-%d %s]",
                spliceJunction.gene(),
                spliceJunction.chromosome(),
                spliceJunction.junctionStart(),
                spliceJunction.junctionEnd(),
                spliceJunction.type());
    }

    // only show transcript ID for non canonical transcripts
    private static String geneTranscriptLabel(String gene, boolean isCanonical, String transcriptId)
    {
        return isCanonical ? gene : String.format("%s(%s)", gene, transcriptId);
    }
}
