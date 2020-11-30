package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class CopyNumberExtractorTest {

    @Test
    public void canExtractCopyNumbersAmps() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HmfGenePanelSupplier.allGenesMap37());

        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC, "any", "AKT1", "AKT1  amp");
        Map<Feature, KnownCopyNumber> featureAmps = Maps.newHashMap();
        featureAmps.put(viccEntry.features().get(0),
                ImmutableKnownCopyNumber.builder().gene("AKT1").type(CopyNumberType.AMPLIFICATION).build());

        assertEquals(featureAmps, copyNumberExtractor.extractAmplificationsDeletions(viccEntry));
    }

    @Test
    public void canExtractCopyNumbersAmpsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HmfGenePanelSupplier.allGenesMap37());

        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC, "any", "TET", "AMPLIFICATION");
        Map<Feature, KnownCopyNumber> featureAmpsUnknown = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertTrue(featureAmpsUnknown.isEmpty());
    }

    @Test
    public void canExtractCopyNumbersDels() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HmfGenePanelSupplier.allGenesMap37());

        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC, "any", "PTEN", "DELETION");
        Map<Feature, KnownCopyNumber> featureDels = Maps.newHashMap();
        featureDels.put(viccEntry.features().get(0), ImmutableKnownCopyNumber.builder().gene("PTEN").type(CopyNumberType.DELETION).build());

        assertEquals(featureDels, copyNumberExtractor.extractAmplificationsDeletions(viccEntry));
    }

    @Test
    public void canExtractCopyNumbersDelsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HmfGenePanelSupplier.allGenesMap37());

        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC, "any", "AC", "DELETION");
        Map<Feature, KnownCopyNumber> featureDelsUnknown = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertTrue(featureDelsUnknown.isEmpty());
    }

}