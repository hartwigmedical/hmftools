package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

public class PeakGnomadFrequenciesCheckerTest extends PurityTestBase
{
    CandidatePeak Level1 = new CandidatePeak(0.1);

    @Before
    public void setup()
    {
        Level1 = new CandidatePeak(0.1);
        Level1.test(evidenceWithDepthAndAltCount(_1, 1000, 100, 90));
        Level1.test(evidenceWithDepthAndAltCount(_1, 2000, 100, 91));
        Level1.test(evidenceWithDepthAndAltCount(_1, 3000, 100, 10));
        Level1.test(evidenceWithDepthAndAltCount(_1, 4000, 100, 9));
    }

    @Test
    public void highAlleleDepthVariantsAreSkewedTest()
    {
        GnomadFrequencySupplier supplier = (chromosome, position) ->
        {
            if(position < 2500)
            {
                return 0.9;
            }
            return 0.5;
        };
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.45));
    }

    @Test
    public void lowAlleleDepthVariantsAreSkewedTest()
    {
        GnomadFrequencySupplier supplier = (chromosome, position) ->
        {
            if(position > 2500)
            {
                return 0.9;
            }
            return 0.5;
        };
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.45));
    }

    @Test
    public void lowAndHighVariantsAreSkewedTest()
    {
        GnomadFrequencySupplier supplier = (chromosome, position) -> 0.9;
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.45));
    }

    @Test
    public void variantMeansMatchExpectedTest()
    {
        GnomadFrequencySupplier supplier = (chromosome, position) -> 0.5;
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertTrue(classifier.checkGnomadFrequencies(supplier, 0.45));
    }

    @Test
    public void checkToleranceTest()
    {
        GnomadFrequencySupplier supplier = (chromosome, position) -> 0.5;
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertTrue(classifier.checkGnomadFrequencies(supplier, 0.36));
        assertTrue(classifier.checkGnomadFrequencies(supplier, 0.64));
    }

    @Test
    public void includeHetAndHomCapturesTest()
    {
        Level1 = new CandidatePeak(0.1);
        Level1.test(evidenceWithDepthAndAltCount(_1, 1000, 1000, 900)); // hom
        Level1.test(evidenceWithDepthAndAltCount(_1, 2000, 1000, 50)); // het
        GnomadFrequencySupplier supplier = (chromosome, position) ->
        {
            if(position < 1500)
            {
                return 0.95;
            }
            return 0.35;
        };
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.95));
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.35));
    }

    @Test
    public void handleNoHighAlleleDepthVariantsTest()
    {
        Level1 = new CandidatePeak(0.1);
        Level1.test(evidenceWithDepthAndAltCount(_1, 2000, 1000, 50)); // het
        GnomadFrequencySupplier supplier = (chromosome, position) -> 0.5;
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.5));
    }

    @Test
    public void handleNoLowAlleleDepthVariantsTest()
    {
        Level1 = new CandidatePeak(0.1);
        Level1.test(evidenceWithDepthAndAltCount(_1, 2000, 1000, 950)); // het
        GnomadFrequencySupplier supplier = (chromosome, position) -> 0.5;
        PeakGnomadFrequenciesChecker classifier = new PeakGnomadFrequenciesChecker(Level1);
        assertFalse(classifier.checkGnomadFrequencies(supplier, 0.5));
    }

}
