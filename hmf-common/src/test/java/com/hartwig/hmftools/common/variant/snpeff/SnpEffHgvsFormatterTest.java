package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SnpEffHgvsFormatterTest {

    @Test
    public void canConvertHgvsCoding() {
        SnpEffAnnotation withCodingAnnotation = SnpEffAnnotationTestFactory.builder().hgvsCoding("c.224C>G").build();
        assertEquals("224C>G", SnpEffHgvsFormatter.formattedHgvsCoding(withCodingAnnotation));

        SnpEffAnnotation withoutCodingAnnotation = SnpEffAnnotationTestFactory.builder().hgvsCoding("hello").build();
        assertEquals("hello", SnpEffHgvsFormatter.formattedHgvsCoding(withoutCodingAnnotation));
    }

    @Test
    public void canConvertHgvsProtein() {
        SnpEffAnnotation snv = SnpEffAnnotationTestFactory.builder().hgvsProtein("p.Val600Lys").build();
        assertEquals("V600K", SnpEffHgvsFormatter.formattedHgvsProtein(snv));

        SnpEffAnnotation frameshift = SnpEffAnnotationTestFactory.builder().hgvsProtein("p.Ser200fs").build();
        assertEquals("S200fs", SnpEffHgvsFormatter.formattedHgvsProtein(frameshift));
    }
}