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
        checkSaavInputResultsInErrorWithMessage("BRAF p.V600E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF:pV600E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF:x.V600E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF:p.V6.00E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF:p.VsixhundredE", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF:p.V600", HGVS_FORMAT_REQUIRED);
        // Note that we allow BRAF:p.600E as the V at 600 is known from the transcript, so is redundant
        checkSaavInputResultsInErrorWithMessage("BRAF p.VV600E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF p.V600EV", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF p.Vali600E", HGVS_FORMAT_REQUIRED);
        checkSaavInputResultsInErrorWithMessage("BRAF p.Val600Isoleucine", HGVS_FORMAT_REQUIRED);

        checkSaavInputResultsInErrorWithMessage("BRUZ:p.A123G", "BRUZ is not a known gene");
        checkSaavInputResultsInErrorWithMessage("BRAF:p.J123G", "Not a valid amino acid identifier: J");
        checkSaavInputResultsInErrorWithMessage("BRAF:p.G123J", "Not a valid amino acid identifier: J");
        checkSaavInputResultsInErrorWithMessage("BRAF:p.G123,456E", HGVS_FORMAT_REQUIRED);
    }

    @Test
    public void canReferToGeneByName()
    {
        SingleAminoAcidVariant variant = variationParser.parseSingleAminoAcidVariant("ZYX:p.Pro46Ala");
        Assert.assertEquals(46, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("ENSG00000159840", variant.Gene.GeneId);
    }

    @Test
    public void canReferToGeneByEnsemblId()
    {
        SingleAminoAcidVariant variant = variationParser.parseSingleAminoAcidVariant("ENSG00000159840:p.Pro46Ala");
        Assert.assertEquals(46, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("ZYX", variant.Gene.GeneName);
    }

    @Test
    public void referenceAminoAcidIsNotRequired()
    {
        SingleAminoAcidVariant variant = variationParser.parseSingleAminoAcidVariant("BRAF:p.600E");
        Assert.assertEquals(600, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void parseSingleAminoAcidVariantTest()
    {
        SingleAminoAcidVariant variant = variationParser.parseSingleAminoAcidVariant("BRAF", "V600E");
        Assert.assertEquals(600, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void aminoAcidNameIsConvertedToSingleLetter()
    {
        SingleAminoAcidVariant variant = variationParser.parseSingleAminoAcidVariant("BRAF:p.Val600Glu");
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void parseDeletionInsertion()
    {
        DeletionInsertion di = variationParser.parseDeletionInsertion("EGFR:p.L747_A750delinsP");
        Assert.assertEquals("EGFR", di.Gene.GeneName);
        Assert.assertEquals(747, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(4, di.changedReferenceSequenceLength());
        Assert.assertEquals("P", di.altAminoAcidSequence());

        // ADCK2 Glu301_Thr303delinsGlnGln, which is E301_T303delinsQQ
        DeletionInsertion di2 = variationParser.parseDeletionInsertion("ADCK2:p.Glu301_Thr303delinsGlnGln");
        Assert.assertEquals("ADCK2", di2.Gene.GeneName);
        Assert.assertEquals(301, di2.positionOfFirstAlteredCodon());
        Assert.assertEquals(3, di2.changedReferenceSequenceLength());
        Assert.assertEquals("QQ", di2.altAminoAcidSequence());
    }

    @Test
    public void parseDeletionInsertionWithGene()
    {
        DeletionInsertion di = variationParser.parseDeletionInsertion("EGFR", "L747_A750delinsP");
        Assert.assertEquals("EGFR", di.Gene.GeneName);
        Assert.assertEquals(747, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(4, di.changedReferenceSequenceLength());
        Assert.assertEquals("P", di.altAminoAcidSequence());
    }

    @Test
    public void deletionInsertionErrors()
    {
        checkDiInputResultsInErrorWithMessage("EGFR:p.L747_A740delinsP","End position must not be before start position");
        checkDiInputResultsInErrorWithMessage("EGFR:p.L747_A750delinsAlaPralineLeu","Not a known amino acid: Praline");
    }

    private void checkSaavInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variationParser.parseSingleAminoAcidVariant(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        Assert.assertEquals(expectedMessage, receivedMessage);
    }

    private void checkDiInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variationParser.parseDeletionInsertion(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        Assert.assertEquals(expectedMessage, receivedMessage);
    }
}
