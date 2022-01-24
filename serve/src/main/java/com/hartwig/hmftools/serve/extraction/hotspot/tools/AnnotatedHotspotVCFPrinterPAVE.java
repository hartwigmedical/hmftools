package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class AnnotatedHotspotVCFPrinterPAVE {

    private static final Logger LOGGER = LogManager.getLogger(AnnotatedHotspotVCFPrinterPAVE.class);

    public static void main(String[] args) throws IOException {
        String annotatedInputVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspots.vcf";
        new AnnotatedHotspotVCFPrinterPAVE().run(annotatedInputVcf);
    }

    public void run(@NotNull String annotatedInputVcf) throws IOException {
        List<HmfTranscriptRegion> canonicalTranscripts = HmfGenePanelSupplier.allGeneList37();

        LOGGER.info("Simplifying variants from '{}'", annotatedInputVcf);
        AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(annotatedInputVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            VariantImpact impact = VariantImpactSerialiser.fromVariantContext(variant);

            String canonicalProtein = impact.CanonicalHgvsProtein;

            if (canonicalProtein.isEmpty()) {
                canonicalProtein = "-";
            }

            StringJoiner joiner = new StringJoiner("|");
            joiner.add(variant.getContig())
                    .add(String.valueOf(variant.getStart()))
                    .add(variant.getAlleles().get(0).getBaseString())
                    .add(variant.getAlleles().get(1).getBaseString())
                    .add(impact.CanonicalGeneName)
                    .add(impact.CanonicalTranscript)
                    .add(impact.CanonicalEffect)
                    .add(impact.CanonicalHgvsCoding)
                    .add(AminoAcids.forceSingleLetterProteinAnnotation(canonicalProtein));

            Object input = variant.getAttribute(VCFWriterFactory.INPUT_FIELD);
            if (input != null) {
                joiner.add(input.toString());
            }

            List<String> sources = variant.getAttributeAsStringList(VCFWriterFactory.SOURCES_FIELD, Strings.EMPTY);
            if (!sources.isEmpty()) {
                StringJoiner sourceJoiner = new StringJoiner(",");
                for (String source : sources) {
                    sourceJoiner.add(source);
                }
                joiner.add(sourceJoiner.toString());
            } else {
                LOGGER.warn("No sources found on {}", variant);
            }

            System.out.println(joiner.toString());
        }
    }
}