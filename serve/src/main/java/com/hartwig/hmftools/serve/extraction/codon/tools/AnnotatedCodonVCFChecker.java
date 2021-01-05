package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.IOException;

import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class AnnotatedCodonVCFChecker {
    private static final Logger LOGGER = LogManager.getLogger(AnnotatedCodonVCFChecker.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedCodonVcf = System.getProperty("user.home") + "/hmf/tmp/annotated_test.vcf";

        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedCodonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            String inputCodon = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|")[2];
            LOGGER.info(SnpEffAnnotationFactory.fromContext(variant).get(0).hgvsProtein());
        }
    }
}
