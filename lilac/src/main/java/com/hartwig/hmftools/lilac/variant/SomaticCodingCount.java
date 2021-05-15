package com.hartwig.hmftools.lilac.variant;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SomaticCodingCount
{
    public final HlaAllele Allele;

    public double InframeIndel;
    public double Missense;
    public double Nonsense;
    public double Splice;
    public double Synonymous;

    public SomaticCodingCount(final HlaAllele allele, double inframeIndel, double missense, double nonsense, double splice, double synonymous)
    {
        Allele = allele;
        InframeIndel = inframeIndel;
        Missense = missense;
        Nonsense = nonsense;
        Splice = splice;
        Synonymous = synonymous;
    }

    public String toString()
    {
        return "SomaticCodingCount(allele=" + Allele + ", inframeIndel=" + InframeIndel + ", missense=" + Missense
                + ", nonsense=" + Nonsense + ", splice=" + Splice + ", synonymous=" + Synonymous + ")";
    }

    public double total() { return InframeIndel + Missense + Nonsense + Splice + Synonymous; }

    public static List<SomaticCodingCount> create(final List<HlaAllele> winners)
    {
        return winners.stream()
                .map(x -> new SomaticCodingCount(x, 0, 0, 0, 0, 0)).collect(Collectors.toList());
    }

    public static void addVariant(
            final List<SomaticCodingCount> codingCounts, final VariantContextDecorator variant, final List<HlaAllele> variantAlleles)
    {
        boolean isIndel = variant.alt().length() != variant.ref().length();
        CodingEffect codingEffect = variant.canonicalCodingEffect();
        addVariant(codingCounts, isIndel, codingEffect, variantAlleles);
    }

    public static void addVariant(
            final List<SomaticCodingCount> codingCounts, boolean isIndel, final CodingEffect effect, final List<HlaAllele> variantAlleles)
    {
        List<SomaticCodingCount> applicableCodingCounts = codingCounts.stream()
                .filter(x -> x.Allele.equals(variantAlleles)).collect(Collectors.toList());

        double contribution = 1.0 / variantAlleles.size();

        for(HlaAllele allele : variantAlleles)
        {
            codingCounts.stream().filter(x -> x.Allele.equals(allele)).forEach(x -> x.addVariant(isIndel, effect, contribution));
        }
        /*
            // Alternative, give only to first to first
            if (counts.isNotEmpty()) {
                result.add(counts[0].addVariant(indel, effect, contribution))
                result.addAll(counts.takeLast(counts.size - 1))
            }
        */
    }

    private void addVariant(boolean indel, CodingEffect effect, double contribution)
    {
        if(indel && effect == CodingEffect.MISSENSE)
        {
            InframeIndel += contribution;
        }

        switch(effect)
        {
            case MISSENSE:
                Missense += contribution;
                break;

            case NONSENSE_OR_FRAMESHIFT:
                Nonsense += contribution;
                break;

            case SPLICE:
                Splice += contribution;
                break;

            case SYNONYMOUS:
                Synonymous += contribution;
                break;

            default:
                break;
        }
    }

}
