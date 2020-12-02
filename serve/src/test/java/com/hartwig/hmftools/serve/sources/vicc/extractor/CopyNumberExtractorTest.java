package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Test;

public class CopyNumberExtractorTest {

    @Test
    public void canExtractCopyNumbersAmps() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("AKT1", "AKT1  amp");

        Map<Feature, KnownCopyNumber> featureAmps = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertEquals(1, featureAmps.size());
        assertEquals("AKT1", featureAmps.get(viccEntry.features().get(0)).gene());
        assertEquals(CopyNumberType.AMPLIFICATION, featureAmps.get(viccEntry.features().get(0)).type());
    }

    @Test
    public void canExtractCopyNumbersAmpsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NOT-A-GENE", "AMPLIFICATION");

        Map<Feature, KnownCopyNumber> featureAmpsUnknown = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertTrue(featureAmpsUnknown.isEmpty());
    }

    @Test
    public void canExtractCopyNumbersDels() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("PTEN", "DELETION");

        Map<Feature, KnownCopyNumber> featureDels = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertEquals(1, featureDels.size());
        assertEquals("PTEN", featureDels.get(viccEntry.features().get(0)).gene());
        assertEquals(CopyNumberType.DELETION, featureDels.get(viccEntry.features().get(0)).type());
    }

    @Test
    public void canExtractCopyNumbersDelsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NOT-A-GENE", "DELETION");

        Map<Feature, KnownCopyNumber> featureDelsUnknown = copyNumberExtractor.extractAmplificationsDeletions(viccEntry);

        assertTrue(featureDelsUnknown.isEmpty());
    }
}