package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.VariationParser.HGVS_FORMAT_REQUIRED;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class VariationParserTest extends TransvalTestBase
{

    private VariationParser variationParser;

    @Before
    public void setUp()
    {
        variationParser = transval.variationParser();
    }

    @Test
    public void errorMessages()
    {
        checkInputResultsInErrorWithMessage("BRAF p.V600E", HGVS_FORMAT_REQUIRED);
        checkInputResultsInErrorWithMessage("BRAF:pV600E", HGVS_FORMAT_REQUIRED);
        checkInputResultsInErrorWithMessage("BRAF:x.V600E", HGVS_FORMAT_REQUIRED);
        checkInputResultsInErrorWithMessage("BRAF:p.V6.00E", HGVS_FORMAT_REQUIRED);
        checkInputResultsInErrorWithMessage("BRAF:p.VsixhundredE", HGVS_FORMAT_REQUIRED);
        checkInputResultsInErrorWithMessage("BRAF:p.V600", HGVS_FORMAT_REQUIRED);
        // Note that we allow BRAF:p.600E as the V at 600 is known from the transcript, so is redundant

        checkInputResultsInErrorWithMessage("BRUZ:p.A123G", "BRUZ is not a known gene");
        checkInputResultsInErrorWithMessage("BRAF:p.J123G", "J does not represent an amino acid");
        checkInputResultsInErrorWithMessage("BRAF:p.G123J", "J does not represent an amino acid");
        checkInputResultsInErrorWithMessage("BRAF:p.G123,456E", HGVS_FORMAT_REQUIRED);
    }

    @Test
    public void canReferToGeneByName()
    {
        SingleAminoAcidVariant variant = variationParser.parse("ZYX:p.Pro46Ala");
        Assert.assertEquals(46, variant.Position);
        Assert.assertEquals("ENSG00000159840", variant.Gene.GeneId);
    }

    @Test
    public void canReferToGeneByEnsemblId()
    {
        SingleAminoAcidVariant variant = variationParser.parse("ENSG00000159840:p.Pro46Ala");
        Assert.assertEquals(46, variant.Position);
        Assert.assertEquals("ZYX", variant.Gene.GeneName);
    }

    @Test
    public void referenceAminoAcidIsNotRequired()
    {
        SingleAminoAcidVariant variant = variationParser.parse("BRAF:p.600E");
        Assert.assertEquals(600, variant.Position);
        Assert.assertEquals("E", variant.Variant);
    }

    @Test
    public void aminoAcidNameIsConvertedToSingleLetter()
    {
        SingleAminoAcidVariant variant = variationParser.parse("BRAF:p.Val600Glu");
        Assert.assertEquals("E", variant.Variant);
    }

    private void checkInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variationParser.parse(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        Assert.assertEquals(expectedMessage, receivedMessage);
    }
}
