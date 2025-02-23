package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.Checks.HGVS_FORMAT_REQUIRED;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class VariantParserTest extends TransvalTestBase
{

    private VariantParser variantParser;

    @Before
    public void setUp()
    {
        variantParser = transval.variationParser();
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
    public void parseDeletionWithGene()
    {
        Deletion di = (Deletion) variantParser.parseVariantForGene("EGFR", "N73_Y74del");
        Assert.assertEquals("EGFR", di.mGene.GeneName);
        Assert.assertEquals(73, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(2, di.mRefLength);
    }

    @Test
    public void parseFrameshiftWithGene()
    {
        Frameshift frameshift = (Frameshift) variantParser.parseVariantForGene("EGFR", "N73fs");
        Assert.assertEquals("EGFR", frameshift.mGene.GeneName);
        Assert.assertEquals(73, frameshift.positionOfFirstAlteredCodon());
        Assert.assertEquals("N", frameshift.mFirstChangedAminoAcid.symbol);
        Assert.assertEquals(1, frameshift.mRefLength);
    }

    @Test
    public void parseStopGainedWithGene()
    {
        StopGained variant = (StopGained) variantParser.parseVariantForGene("EGFR", "N73*");
        Assert.assertEquals("EGFR", variant.mGene.GeneName);
        Assert.assertEquals(73, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("N", variant.mFirstChangedAminoAcid.symbol);
        Assert.assertEquals(1, variant.mRefLength);
    }

    @Test
    public void parseInsertionWithGene()
    {
        Insertion di = (Insertion) variantParser.parseInsertion("EGFR", "N73_Y74insSPQR");
        Assert.assertEquals("EGFR", di.mGene.GeneName);
        Assert.assertEquals(74, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(2, di.mRefLength);
        Assert.assertEquals("SPQR", di.mInsertedSequence.sequence());
    }

    @Test
    public void canReferToGeneByName()
    {
        SingleAminoAcidVariant variant = variantParser.parseSingleAminoAcidVariant("ZYX:p.Pro46Ala");
        Assert.assertEquals(46, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("ENSG00000159840", variant.mGene.GeneId);
    }

    @Test
    public void canReferToGeneByEnsemblId()
    {
        SingleAminoAcidVariant variant = variantParser.parseSingleAminoAcidVariant("ENSG00000159840:p.Pro46Ala");
        Assert.assertEquals(46, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("ZYX", variant.mGene.GeneName);
    }

    @Test
    public void referenceAminoAcidIsNotRequired()
    {
        SingleAminoAcidVariant variant = variantParser.parseSingleAminoAcidVariant("BRAF:p.600E");
        Assert.assertEquals(600, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void parseSingleAminoAcidVariantTest()
    {
        SingleAminoAcidVariant variant = variantParser.parseSingleAminoAcidVariant("BRAF", "V600E");
        Assert.assertEquals(600, variant.positionOfFirstAlteredCodon());
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void aminoAcidNameIsConvertedToSingleLetter()
    {
        SingleAminoAcidVariant variant = variantParser.parseSingleAminoAcidVariant("BRAF:p.Val600Glu");
        Assert.assertEquals("E", variant.altValue());
    }

    @Test
    public void parseVariantForGeneTest()
    {
        ProteinVariant pv1 = variantParser.parseVariantForGene("BRAF", "Val600Glu");
        Assert.assertEquals(600, pv1.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv1 instanceof SingleAminoAcidVariant);

        ProteinVariant pv2 = variantParser.parseVariantForGene("ADCK2", "Glu301_Thr303delinsGlnGln");
        Assert.assertEquals(301, pv2.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv2 instanceof DeletionInsertion);

        ProteinVariant pv3 = variantParser.parseVariantForGene("EGFR", "L747_A750del");
        Assert.assertEquals(747, pv3.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv3 instanceof Deletion);

        ProteinVariant pv4 = variantParser.parseVariantForGene("PIK3R1", "Y452dup");
        Assert.assertEquals(452, pv4.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv4 instanceof Duplication);
    }

    @Test
    public void parseDeletionInsertion()
    {
        DeletionInsertion di = (DeletionInsertion) variantParser.parseDeletionInsertion("EGFR:p.L747_A750delinsP");
        Assert.assertEquals("EGFR", di.mGene.GeneName);
        Assert.assertEquals(747, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(4, di.mRefLength);
        Assert.assertEquals("P", di.altAminoAcidSequence());

        // ADCK2 Glu301_Thr303delinsGlnGln, which is E301_T303delinsQQ
        DeletionInsertion di2 = (DeletionInsertion) variantParser.parseDeletionInsertion("ADCK2:p.Glu301_Thr303delinsGlnGln");
        Assert.assertEquals("ADCK2", di2.mGene.GeneName);
        Assert.assertEquals(301, di2.positionOfFirstAlteredCodon());
        Assert.assertEquals(3, di2.mRefLength);
        Assert.assertEquals("QQ", di2.altAminoAcidSequence());
    }

    @Test
    public void parseDeletionInsertionWithGene()
    {
        DeletionInsertion di = (DeletionInsertion) variantParser.parseDeletionInsertion("EGFR", "L747_A750delinsP");
        Assert.assertEquals("EGFR", di.mGene.GeneName);
        Assert.assertEquals(747, di.positionOfFirstAlteredCodon());
        Assert.assertEquals(4, di.mRefLength);
        Assert.assertEquals("P", di.altAminoAcidSequence());
    }

    @Test
    public void deletionInsertionErrors()
    {
        checkDiInputResultsInErrorWithMessage("EGFR:p.L747_A740delinsP", "End position must not be before start position");
        checkDiInputResultsInErrorWithMessage("EGFR:p.L747_A750delinsAlaPralineLeu", "Not a valid amino acid identifier: Praline");
    }

    @Test
    public void parseDuplicationWithGene()
    {
        Duplication dup = (Duplication) variantParser.parseDuplication("PIK3R1", "Y452dup");
        Assert.assertEquals("PIK3R1", dup.mGene.GeneName);
        Assert.assertEquals(452, dup.positionOfFirstAlteredCodon());
        Assert.assertEquals(1, dup.mRefLength);

        Duplication dup2 = (Duplication) variantParser.parseDuplication("PIK3R1", "E458_Y463dup");
        Assert.assertEquals("PIK3R1", dup2.mGene.GeneName);
        Assert.assertEquals(458, dup2.positionOfFirstAlteredCodon());
        Assert.assertEquals(6, dup2.mRefLength);
    }

    private void checkSaavInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variantParser.parseSingleAminoAcidVariant(input);
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
            variantParser.parseDeletionInsertion(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        Assert.assertEquals(expectedMessage, receivedMessage);
    }
}
