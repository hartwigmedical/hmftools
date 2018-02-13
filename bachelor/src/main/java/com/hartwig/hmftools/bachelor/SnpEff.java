package com.hartwig.hmftools.bachelor;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

// this class is logically keyed by TranscriptID (aka IsoForm) and Allele
// there can be many entruies per variant (VCF line) annotated by SnpEff
@Deprecated
public class SnpEff {

    private final String GeneName;
    private final String Transcript;
    private final List<String> Effects;
    private final String AllEffects;
    private final String HGVSc;
    private final String HGVSp;
    private final List<Integer> ProteinPosition;

    private final String allele;

    private SnpEff(final String name, final String allele, final String transcript, final String allEffects, final String hgvsc, final String hgvsp,
            final List<Integer> proteinPosition) {
        GeneName = name;
        this.allele = allele;
        Transcript = transcript;
        Effects = Arrays.asList(allEffects.split("\\&"));
        AllEffects = allEffects;
        HGVSc = hgvsc;
        HGVSp = hgvsp;
        ProteinPosition = proteinPosition;
    }

    // Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
    // eg C|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000351666|protein_coding|11/19|c.1288A>G|p.Ser430Gly|1307/3444|1288/2043|430/680
    private enum Field {
        ALLELE(0),
        ANNOTATION(1),
        ANNOTATION_IMPACT(2),
        GENE_NAME(3),
        GENE_ID(4),
        FEATURE_TYPE(5),
        FEATURE_ID(6),
        TRANSCRIPT_BIOTYPE(7),
        RANK(8),
        HGVSc(9),
        HGVSp(10),
        cDNA(11),
        CDS(12),
        PROTEIN(13),
        MAX(14);

        private final int value;

        Field(final int v) {
            value = v;
        }
    }

    @Nullable
    @Deprecated
    static SnpEff parseAnnotation(final List<String> annotation) {
        if (annotation.size() < Field.MAX.value) {
            return null;
        }
        if (annotation.get(Field.FEATURE_TYPE.value).equals("transcript")) {
            final String geneName = annotation.get(Field.GENE_NAME.value);
            String allele = annotation.get(Field.ALLELE.value).trim();

            if(allele.charAt(0) == '[') {
                allele = allele.substring(1);
            }

            final String transcript = annotation.get(Field.FEATURE_ID.value);
            //final List<String> effects = Arrays.asList(annotation.get(Field.ANNOTATION.value).split("\\&"));
            final String allEffects = annotation.get(Field.ANNOTATION.value);
            final String hgvsc = annotation.get(Field.HGVSc.value).replaceFirst("^c\\.", "");
            final String hgvsp = annotation.get(Field.HGVSp.value).replaceFirst("^p\\.", "");

            final List<Integer> proteinPosition = Arrays.stream(annotation.get(Field.PROTEIN.value).split("/"))
                    .filter(s -> !s.isEmpty())
                    .map(Integer::parseInt)
                    .collect(Collectors.toList());

            return new SnpEff(geneName, allele, transcript, allEffects, hgvsc, hgvsp, proteinPosition);
        }
        return null;
    }

    public String allele() {
        return this.allele;
    }

    public String gene() {
        return GeneName;
    }

    public String transcript() {
        return Transcript;
    }

    public List<String> effects() {
        return Effects;
    }

    public String allEffects() {
        return AllEffects;
    }

    public String hgvsCoding() {
        return HGVSc;
    }

    public String hgvsProtein() {
        return HGVSp;
    }

    public List<Integer> proteinPosition() {
        return ProteinPosition;
    }
}
