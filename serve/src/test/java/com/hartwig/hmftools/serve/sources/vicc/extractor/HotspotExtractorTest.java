package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.transvar.Transvar;
import com.hartwig.hmftools.serve.transvar.TransvarTest;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class HotspotExtractorTest {

    @Test
    public void canExtractHotspot() {

        HotspotExtractor hotspotExtractor = new HotspotExtractor(ProteinResolverFactory.dummy(), new ProteinAnnotationExtractor());
        Transvar transvar = TransvarTest.returnsSingleTransvarRecord(TransvarTest.createTestRecord());

        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB, "ENST00000261584", "PALB2", "L939W");
        Feature feature = viccEntry.features().get(0);
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();

        hotspotsPerFeature.put(feature,
                transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                        viccEntry.transcriptId(),
                        hotspotExtractor.extractProteinAnnotation(feature)));

    }
}