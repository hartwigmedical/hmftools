package com.hartwig.hmftools.bachelor;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

public class SnpEff {

    final String GeneName;
    final String Transcript;
    final List<String> Effects;
    final String HGVSc;
    final String HGVSp;
    final List<Integer> ProteinPosition;

    private SnpEff(final String name, final String transcript, final List<String> effects, final String hgvsc, final String hgvsp,
            final List<Integer> proteinPosition) {
        GeneName = name;
        Transcript = transcript;
        Effects = effects;
        HGVSc = hgvsc;
        HGVSp = hgvsp;
        ProteinPosition = proteinPosition;
    }

    // Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
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
    static SnpEff parseAnnotation(final List<String> annotation) {
        if (annotation.size() < Field.MAX.value + 1) {
            return null;
        }
        if (annotation.get(Field.FEATURE_TYPE.value).equals("transcript")) {
            final String geneName = annotation.get(Field.GENE_NAME.value);
            final String transcript = annotation.get(Field.FEATURE_ID.value);
            final List<String> effects = Arrays.asList(annotation.get(Field.ANNOTATION.value).split("\\&"));
            final String hgvsc = annotation.get(Field.HGVSc.value).replaceFirst("^c\\.", "");
            final String hgvsp = annotation.get(Field.HGVSp.value).replaceFirst("^p\\.", "");

            final List<Integer> proteinPosition = Arrays.stream(annotation.get(Field.PROTEIN.value).split("/"))
                    .filter(s -> !s.isEmpty())
                    .map(Integer::parseInt)
                    .collect(Collectors.toList());

            return new SnpEff(geneName, transcript, effects, hgvsc, hgvsp, proteinPosition);
        }
        return null;
    }
}
