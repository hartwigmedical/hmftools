package com.hartwig.hmftools.pavereverse.parse;

import static com.hartwig.hmftools.pavereverse.util.Checks.HGVS_FORMAT_REQUIRED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.protein.Deletion;
import com.hartwig.hmftools.pavereverse.protein.DeletionInsertion;
import com.hartwig.hmftools.pavereverse.protein.Duplication;
import com.hartwig.hmftools.pavereverse.protein.Frameshift;
import com.hartwig.hmftools.pavereverse.protein.Insertion;
import com.hartwig.hmftools.pavereverse.protein.ProteinVariant;
import com.hartwig.hmftools.pavereverse.protein.SingleAminoAcidVariant;
import com.hartwig.hmftools.pavereverse.protein.StartLost;
import com.hartwig.hmftools.pavereverse.protein.StopGained;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class ProteinVariantParserTest extends ReversePaveTestBase
{

    private ProteinVariantParser variantParser;

    @Before
    public void setUp()
    {
        variantParser = reversePave.proteinVariantParser();
    }

    @Test
    public void returnVariantForSuppliedTranscriptId()
    {
        ProteinVariant variant = variantParser.parseGeneVariant("ZYX", "ENST00000436448","S2V");
        assertEquals("ZYX", variant.Gene.GeneName);
        assertEquals("ENST00000436448", variant.Transcript.TransName);
    }

    @Test
    public void reportErrorIfNoTranscriptFoundForId()
    {
        String receivedMessage = null;
        try
        {
            variantParser.parseGeneVariant("ZYX", "Whatever", "F2V");
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals("No transcript found. Gene: ENSG00000159840, transcript id: Whatever", receivedMessage);
    }

    @Test
    public void reportErrorIfIdentifiedTranscriptDoesNotMatchSpecification()
    {
        String receivedMessage = null;
        try
        {
            variantParser.parseGeneVariant("ZYX", "ENST00000436448","W2V");
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals("Transcript does not match. Gene: ENSG00000159840, transcript: ENST00000436448, ref: [2,W]_[2,W]", receivedMessage);
    }

    @Test
    public void returnMultipleVariantsIfCanonicalDoesNotMatchSpecification()
    {
        Set<ProteinVariant> variants = variantParser.parseGeneVariants("ZYX", "S2V");
        // See ensembl_trans_amino_acids.csv
        assertEquals(3, variants.size());
        Set<String> transcriptNames = variants.stream().map(variant -> variant.Transcript.TransName).collect(Collectors.toSet());
        assertTrue(transcriptNames.contains("ENST00000436448"));
        assertTrue(transcriptNames.contains("ENST00000446634"));
        assertTrue(transcriptNames.contains("ENST00000449630"));
    }

    @Test
    public void returnJustTheCanonicalTranscriptIfItMatchesTheSpecification()
    {
        Set<ProteinVariant> variants = variantParser.parseGeneVariants("ZYX", "A2V");
        // See ensembl_trans_amino_acids.csv
        assertEquals(1, variants.size());
        ProteinVariant variant = variants.iterator().next();
        assertEquals("ENST00000322764", variant.Transcript.TransName);
    }

    @Test
    public void reportErrorIfNoMatchingTranscriptsFound()
    {
        String receivedMessage = null;
        try
        {
            variantParser.parseGeneVariants("ZYX", "F2V");
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals("No transcript found for gene: ENSG00000159840 matching ref: [2,F]_[2,F]", receivedMessage);
    }

    @Test
    public void reportErrorIfSingleTranscriptExpectedButMultipleNonCanonicalTranscriptsFound()
    {
        String receivedMessage = null;
        try
        {
            variantParser.parseGeneVariant("ZYX", "S2V");
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals("No canonical transcript, but multiple non-canonical transcripts, found for gene: ENSG00000159840, ref: [2,S]_[2,S]", receivedMessage);
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
        Deletion di = (Deletion) variantParser.parseGeneVariant("EGFR", "N73_Y74del");
        assertEquals("EGFR", di.Gene.GeneName);
        assertEquals(73, di.positionOfFirstAlteredCodon());
        assertEquals(2, di.RefLength);
    }

    @Test
    public void parseFrameshiftWithGene()
    {
        Frameshift frameshift = (Frameshift) variantParser.parseGeneVariant("EGFR", "N73fs");
        assertEquals("EGFR", frameshift.Gene.GeneName);
        assertEquals(73, frameshift.positionOfFirstAlteredCodon());
        assertEquals("N", frameshift.FirstChangedAminoAcid.Symbol);
        assertEquals(1, frameshift.RefLength);
    }

    @Test
    public void parseStopGainedWithGene()
    {
        StopGained variant = (StopGained) variantParser.parseGeneVariant("EGFR", "N73*");
        assertEquals("EGFR", variant.Gene.GeneName);
        assertEquals(73, variant.positionOfFirstAlteredCodon());
        assertEquals("N", variant.FirstChangedAminoAcid.Symbol);
        assertEquals(1, variant.RefLength);
    }

    @Test
    public void parseStartLostWithGene()
    {
        StartLost variant = (StartLost) variantParser.parseGeneVariant("KIT", "M1?");
        assertEquals("KIT", variant.Gene.GeneName);
        assertEquals(1, variant.positionOfFirstAlteredCodon());
    }

    @Test
    public void parseInsertionWithGene()
    {
        Insertion di = (Insertion) variantParser.parseGeneVariant("EGFR", "N73_Y74insSPQR");
        assertEquals("EGFR", di.Gene.GeneName);
        assertEquals(74, di.positionOfFirstAlteredCodon());
        assertEquals(2, di.RefLength);
        assertEquals("SPQR", di.InsertedSequence.sequence());
    }

    @Test
    public void canReferToGeneByName()
    {
        SingleAminoAcidVariant variant = (SingleAminoAcidVariant) variantParser.parse("ZYX:p.Pro46Ala");
        assertEquals(46, variant.positionOfFirstAlteredCodon());
        assertEquals("ENSG00000159840", variant.Gene.GeneId);
    }

    @Test
    public void canReferToGeneByEnsemblId()
    {
        SingleAminoAcidVariant variant = (SingleAminoAcidVariant) variantParser.parse("ENSG00000159840:p.Pro46Ala");
        assertEquals(46, variant.positionOfFirstAlteredCodon());
        assertEquals("ZYX", variant.Gene.GeneName);
    }

    @Test
    public void parseSingleAminoAcidVariantTest()
    {
        SingleAminoAcidVariant variant = (SingleAminoAcidVariant) variantParser.parseGeneVariant("BRAF", "V600E");
        assertEquals(600, variant.positionOfFirstAlteredCodon());
        assertEquals("E", variant.altValue());
    }

    @Test
    public void aminoAcidNameIsConvertedToSingleLetter()
    {
        SingleAminoAcidVariant variant = (SingleAminoAcidVariant) variantParser.parse("BRAF:p.Val600Glu");
        assertEquals("E", variant.altValue());
    }

    @Test
    public void parseGeneVariantTest()
    {
        ProteinVariant pv1 = variantParser.parseGeneVariant("BRAF", "Val600Glu");
        assertEquals(600, pv1.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv1 instanceof SingleAminoAcidVariant);

        ProteinVariant pv2 = variantParser.parseGeneVariant("ADCK2", "Glu301_Thr303delinsGlnGln");
        assertEquals(301, pv2.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv2 instanceof DeletionInsertion);

        ProteinVariant pv3 = variantParser.parseGeneVariant("EGFR", "L747_A750del");
        assertEquals(747, pv3.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv3 instanceof Deletion);

        ProteinVariant pv4 = variantParser.parseGeneVariant("PIK3R1", "Y452dup");
        assertEquals(452, pv4.positionOfFirstAlteredCodon());
        Assert.assertTrue(pv4 instanceof Duplication);
    }

    @Test
    public void parseDeletionInsertion()
    {
        DeletionInsertion di = (DeletionInsertion) variantParser.parse("EGFR:p.L747_A750delinsP");
        assertEquals("EGFR", di.Gene.GeneName);
        assertEquals(747, di.positionOfFirstAlteredCodon());
        assertEquals(4, di.RefLength);
        assertEquals(aaSeq("P"), di.Alt);

        // ADCK2 Glu301_Thr303delinsGlnGln, which is E301_T303delinsQQ
        DeletionInsertion di2 = (DeletionInsertion) variantParser.parse("ADCK2:p.Glu301_Thr303delinsGlnGln");
        assertEquals("ADCK2", di2.Gene.GeneName);
        assertEquals(301, di2.positionOfFirstAlteredCodon());
        assertEquals(3, di2.RefLength);
        assertEquals(aaSeq("QQ"), di2.Alt);
    }

    @Test
    public void parseDeletionInsertionWithGene()
    {
        DeletionInsertion di = (DeletionInsertion) variantParser.parseGeneVariant("EGFR", "L747_A750delinsP");
        assertEquals("EGFR", di.Gene.GeneName);
        assertEquals(747, di.positionOfFirstAlteredCodon());
        assertEquals(4, di.RefLength);
        assertEquals(aaSeq("P"), di.Alt);
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
        Duplication dup = (Duplication) variantParser.parseGeneVariant("PIK3R1", "Y452dup");
        assertEquals("PIK3R1", dup.Gene.GeneName);
        assertEquals(452, dup.positionOfFirstAlteredCodon());
        assertEquals(1, dup.RefLength);

        Duplication dup2 = (Duplication) variantParser.parseGeneVariant("PIK3R1", "E458_Y463dup");
        assertEquals("PIK3R1", dup2.Gene.GeneName);
        assertEquals(458, dup2.positionOfFirstAlteredCodon());
        assertEquals(6, dup2.RefLength);
    }

    private void checkSaavInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variantParser.parse(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals(expectedMessage, receivedMessage);
    }

    private void checkDiInputResultsInErrorWithMessage(String input, String expectedMessage)
    {
        String receivedMessage = null;
        try
        {
            variantParser.parse(input);
        }
        catch(final IllegalArgumentException e)
        {
            receivedMessage = e.getMessage();
        }
        assertEquals(expectedMessage, receivedMessage);
    }
}
