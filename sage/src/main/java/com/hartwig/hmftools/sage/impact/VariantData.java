package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.variant.VariantConsequence.NON_CODING_TRANSCRIPT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.consequencesToString;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.VariantType.MNP;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.DELIM;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.ITEM_DELIM;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    private final List<VariantTransImpact> mImpacts;

    private final List<SnpEffAnnotation> mSnpEffAnnotations;

    public VariantData(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        mImpacts = Lists.newArrayList();
        mSnpEffAnnotations = Lists.newArrayList();
    }

    public static VariantData fromContext(final VariantContext variantContext)
    {
        int variantPosition = (int)variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = variantContext.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));

        return new VariantData(chromosome, variantPosition, ref, alt);
    }

    public VariantType type()
    {
        if(Ref.equals(Alt))
            return Ref.length() == 1 ? SNP : MNP;

        return INDEL;
    }

    public List<VariantTransImpact> getImpacts() { return mImpacts; }

    public void addImpact(final VariantTransImpact impact)
    {
        if(impact.TransData == null)
        {
            if(mImpacts.stream().anyMatch(x -> x.Consequence == impact.Consequence))
                return;
        }
        else
        {
            if(mImpacts.stream().filter(x -> x.TransData != null) .anyMatch(x -> x.TransData.TransId == impact.TransData.TransId))
                return;
        }

        mImpacts.add(impact);
    }

    public VariantTransImpact getWorstImpact()
    {
        return null;
    }

    public VariantTransImpact getCanonicalImpact()
    {
        return mImpacts.stream().filter(x -> x.TransData.IsCanonical).findFirst().orElse(null);
    }

    public void setSnpEffAnnotations(final List<SnpEffAnnotation> annotations) { mSnpEffAnnotations.addAll(annotations); };

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s>%s)", Chromosome, Position, Ref, Alt);
    }

    // output methods
    public static String csvHeader()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(csvCommonHeader());
        // sj.add("Chromosome").add("Position").add("Type").add("Ref").add("Alt");
        return sj.toString();
    }

    public String csvData(final EnsemblGeneData geneData)
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(csvCommonData(geneData));
        return sj.toString();
    }

    private static String csvCommonHeader()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add("Chromosome").add("Position").add("Type").add("Ref").add("Alt");
        sj.add("GeneId").add("GeneName");
        return sj.toString();
    }

    private String csvCommonData(final EnsemblGeneData geneData)
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(Chromosome);
        sj.add(String.valueOf(Position));
        sj.add(String.valueOf(type()));
        sj.add(Ref);
        sj.add(Alt);

        if(geneData != null)
        {
            sj.add(geneData.GeneId);
            sj.add(geneData.GeneName);
        }
        else
        {
            sj.add("NO_GENE");
            sj.add("NO_GENE");
        }

        return sj.toString();
    }

    public static String csvTranscriptHeader()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(csvCommonHeader());
        sj.add("TransId");
        sj.add("Consequence");
        sj.add("SeConsequence");
        return sj.toString();
    }

    public List<String> csvTranscriptData(final EnsemblGeneData geneData)
    {
        List<String> transcriptLines = Lists.newArrayList();

        List<SnpEffAnnotation> matchedAnnotations = Lists.newArrayList();

        for(VariantTransImpact impact : mImpacts)
        {
            if(impact.TransData == null)
                continue;

            SnpEffAnnotation annotation = impact.findMatchingAnnotation(mSnpEffAnnotations);

            StringJoiner sj = new StringJoiner(DELIM);
            sj.add(csvCommonData(geneData));

            sj.add(impact.TransData.TransName);
            sj.add(String.valueOf(impact.Consequence));

            if(annotation != null)
            {
                sj.add(consequencesToString(annotation.consequences(), ITEM_DELIM));
                matchedAnnotations.add(annotation);
            }
            else
            {
                sj.add("UNMATCHED");
            }

            transcriptLines.add(sj.toString());
        }

        for(SnpEffAnnotation annotation : mSnpEffAnnotations)
        {
            if(matchedAnnotations.contains(annotation))
                continue;

            if(annotation.consequences().contains(NON_CODING_TRANSCRIPT_VARIANT))
                continue;

            StringJoiner sj = new StringJoiner(DELIM);
            sj.add(csvCommonData(geneData));

            sj.add(annotation.featureID());
            sj.add("UNMATCHED");
            sj.add(consequencesToString(annotation.consequences(), ITEM_DELIM));
            transcriptLines.add(sj.toString());
        }

        return transcriptLines;
    }

}
