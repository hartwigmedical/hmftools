package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.transvar.Transvar;
import com.hartwig.hmftools.serve.transvar.TransvarTest;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Ignore;
import org.junit.Test;

public class HotspotExtractorTest {

    @Test
    @Ignore
    public void canExtractHotspot() {
        HotspotExtractor hotspotExtractor = new HotspotExtractor(ProteinResolverFactory.dummy(), new ProteinAnnotationExtractor());

        ViccEntry viccEntry =
                ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB, "ENST00000261584", "PALB2", "L939W", "7", "10", null);
        Feature feature = viccEntry.features().get(0);
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();
        Transvar transvar = TransvarTest.returnsSingleTransvarRecord(TransvarTest.createTestRecordHotspot(feature.chromosome(),
                Integer.valueOf(feature.start()),
                viccEntry.transcriptId()));

        hotspotsPerFeature.put(feature,
                transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                        viccEntry.transcriptId(),
                        hotspotExtractor.extractProteinAnnotation(feature)));

        Map<Feature, List<VariantHotspot>> hotspotsPerFeatureExpect = Maps.newHashMap();
        hotspotsPerFeatureExpect.put(ImmutableFeature.builder()
                .name("L939W")
                .chromosome("7")
                .start("10")
                .provenance(Lists.newArrayList())
                .geneSymbol("PALB2")
                .synonyms(Lists.newArrayList())
                .links(Lists.newArrayList())
                .build(), Lists.newArrayList());

        assertEquals(hotspotsPerFeatureExpect, hotspotsPerFeature);

    }
}