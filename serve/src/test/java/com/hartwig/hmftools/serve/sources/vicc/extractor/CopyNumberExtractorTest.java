package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneCheckerTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Test;

public class CopyNumberExtractorTest {

    private static final GeneChecker HG19_GENE_CHECKER = GeneCheckerTestFactory.buildForHG19();

    @Test
    public void canExtractCopyNumbersAmps() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("AKT1", "AKT1  amp");

        Map<Feature, KnownCopyNumber> featureAmps = copyNumberExtractor.extract(viccEntry);

        assertEquals(1, featureAmps.size());
        assertEquals("AKT1", featureAmps.get(viccEntry.features().get(0)).gene());
        assertEquals(CopyNumberType.AMPLIFICATION, featureAmps.get(viccEntry.features().get(0)).type());
    }

    @Test
    public void canExtractCopyNumbersAmpsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NOT-A-GENE", "AMPLIFICATION");

        Map<Feature, KnownCopyNumber> featureAmpsUnknown = copyNumberExtractor.extract(viccEntry);

        assertTrue(featureAmpsUnknown.isEmpty());
    }

    @Test
    public void canExtractCopyNumbersDels() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("PTEN", "DELETION");

        Map<Feature, KnownCopyNumber> featureDels = copyNumberExtractor.extract(viccEntry);

        assertEquals(1, featureDels.size());
        assertEquals("PTEN", featureDels.get(viccEntry.features().get(0)).gene());
        assertEquals(CopyNumberType.DELETION, featureDels.get(viccEntry.features().get(0)).type());
    }

    @Test
    public void canExtractCopyNumbersDelsUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NOT-A-GENE", "DELETION");

        Map<Feature, KnownCopyNumber> featureDelsUnknown = copyNumberExtractor.extract(viccEntry);

        assertTrue(featureDelsUnknown.isEmpty());
    }
}